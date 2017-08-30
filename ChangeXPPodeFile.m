function success = ChangeXPPodeFile(filename,parset)
%ChangeXPPodeFile  Change parameters / initial conditions in XPP ODE file
%  (c) Robert Clewley, August 2004
%
%  Usage: success = ChangeXPPodeFile(filename,parset) changes the parameters and initial
%  conditions named in the parset structure to new values by directly modifying the named
%  ODE file. The return value success is a boolean.
%
%  IMPORTANT: All initial conditions to be changed must be explicitly declared in the ODE file
%  using 'init' or 'i', not using the format 'x(0)' for a variable x, etc. Also, implicitly-declared
%  zero initial conditions cannot be altered by this function. For these, you might find it better
%  to use the ChangeXPPsetFile function instead, with the setfile option of RunXPP.
%
%  'number' declarations in XPP are also not supported as range parameters.
%
%  Input arguments:
%  filename should end with the extension '.ode'
%  parset is a structure array with fields 'type', 'name', and 'val'. type must be one
%   of the strings 'PAR' or 'IC'.

success = false;

if isempty(filename) | exist(filename,'file') ~= 2
    disp(' Problem with filename passed: file does not exist in specified path')
    return
else
    if ~strcmp('.ode',filename(length(filename)-3:length(filename)))
        disp(' Problem with filename passed: must end in `.ode`')
        return
    end
end

fid = fopen(filename,'r');
file_done = false;
pars_found = false;
lineCount = 0;
odefile = {};
lenpars = length(parset);
parnames = cell(1,lenpars);
parnames = {parset(1:lenpars).name}; % make a simple cell array list of parnames to be changed
foundOccurrence = zeros(1,lenpars);

fieldsexpected = {'name','val','type'};
for ix=1:lenpars
    if ~isempty(setdiff(fieldnames(parset(ix)),fieldsexpected))
        fprintf(' Incorrect fields found in argument parset for entry %d\n',ix)
        fclose(fid);
        return
    end
	if ~ismember(parset(ix).type,{'PAR','IC'})
        fprintf(' Invalid type, %s, specified in entry %d\n',parset(ix).type,ix)
        fclose(fid);
        return
	end
end

% change specified pars when encountered
parsDone = false;
while ~parsDone
    fline = fgetl(fid);
    if ~ischar(fline) % then premature EOF
        disp(' End of file reached before specified params found. Cannot continue')
        fclose(fid);
        return
    end
    [token rest] = strtok(fline);
    isParType = strcmp(token, 'par') | strcmp(token, 'p');
    isICType = strcmp(token, 'init') | strcmp(token, 'i');
    if isParType | isICType
        [parsed numparsed] = ParseParLine(rest);
        if numparsed == 0 || parsed(numparsed).num == NaN
            disp(' Problem parsing numerical value of parameter / initial condition in ODE file. Cannot continue')
            fclose(fid);
            return
        end
        if isParType
            flineNew = ['par '];
            typeStr = 'par';
        else
            flineNew = ['init '];
            typeStr = 'ic';
        end
        for i=1:numparsed
            [ispres ix] = ismember(parsed(i).name, parnames);
            if ispres % then a par to be changed has been found
                if (strcmp(parset(ix).type,'PAR') & ~isParType) | (strcmp(parset(ix).type,'IC') & isParType)
                    fprintf(' Found %s `%s` that didn`t correspond to specified type %s\n',typeStr,parset(ix).name,parset(ix).type)
                    fclose(fid);
                    return
                end
                foundOccurrence(ix) = 1;
                if numparsed==1 | i==numparsed
                    flineNew = [ flineNew, parsed(i).name '=' num2str(parset(ix).val)  ];
                else
                    flineNew = [ flineNew, parsed(i).name '=' num2str(parset(ix).val) ', ' ];
                end
            else
                if numparsed==1 | i==numparsed
                    flineNew = [ flineNew, parsed(i).name '=' num2str(parsed(i).num)  ];
                else
                    flineNew = [ flineNew, parsed(i).name '=' num2str(parsed(i).num) ',' ];
                end
            end
        end
    else
        flineNew = fline;
    end
    odefile = {odefile{:}, flineNew};
    lineCount = lineCount + 1;
    if sum(foundOccurrence) == lenpars % then made all changes necessary
        parsDone = true;
    end
end

% pull in rest of file once params changed
while ~file_done
    fline = fgetl(fid);
    if ~ischar(fline) % then EOF
        file_done = true;
    else
        odefile = {odefile{:}, fline};
        lineCount = lineCount + 1;
    end
end
fclose(fid);

fidw = fopen( filename ,'w');
if fidw ~= -1
	for line=1:lineCount
        fprintf(fidw, '%s\r\n', odefile{line});
	end
	fclose(fidw);
else
    fprintf(' Problem opening file %s for writing. Cannot continue\n',filename);
    return
end
success = true;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [parsed, numparsed] = ParseParLine(full_line)
parseddef.num = NaN;
parseddef.name = '';
fline = full_line(isspace(full_line)==0); % get rid of whitespace
fline = strrep(fline,'=',' '); % convert = to space
fline = strrep(fline,',',' '); % convert , to space
numparsed = 0;
endofline = false;
while ~endofline
    [nameStr restStr] = strtok(fline);
    if isempty(nameStr)
        endofline = true;
        continue
    end
    numparsed = numparsed + 1;
    [numStr restStr] = strtok(restStr);
    if isempty(numStr) | isempty(nameStr)
        disp('Param parse error: Parameter name or value missing!')
        parsed = parseddef;
        return
    end
    if ~isNum(numStr(isspace(numStr)==0))
        disp('Param parse error: Parameter value not a number!')
        parsed = parseddef;
        return
    end
    parsed(numparsed).name = nameStr(isspace(nameStr)==0); % get rid of all whitespace
    parsed(numparsed).num  = str2num(numStr);
    fline = restStr;
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = isNum(data)
% parameter `data` is a string
result = false;
pointCount = 0;
num = [char(48:57),'.','-','e'];
for i=1:length(data)
    if data(i) == '.'
        pointCount = pointCount + 1;
        if pointCount > 1
            return % false
        end
    end
    if ~ismember(data(i),num)
        return % false
    end
end

result = true; % can only get here if all chars were numeric
return
