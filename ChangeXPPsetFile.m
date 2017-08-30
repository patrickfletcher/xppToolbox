function success = ChangeXPPsetFile(filename,parset)
%ChangeXPPsetFile  Change parameters / initial conditions in XPP SET file
%  (c) Robert Clewley, August 2004
%
%  Usage: success = ChangeXPPsetFile(filename,parset) changes the values of the
%  parameters and initial conditions named in the parset structure to by directly
%  modifying the named SET file. The return value success is a boolean.
%
%  Input arguments:
%  filename should end with the extension '.set'
%  parset is a structure array with fields 'type', 'name', and 'val'. type must be one
%   of the strings 'PAR' or 'IC'.


success = false;

if isempty(filename) | exist(filename,'file') ~= 2
    disp(' Problem with filename passed: file does not exist in specified path')
    return
else
    if ~strcmp('.set',filename(length(filename)-3:length(filename)))
        disp(' Problem with filename passed: must end in `.set`')
        return
    end
end

fid = fopen(filename,'r');
lenpars = length(parset);
file_done = false;
ics_found = false;
pars_found = false;
lineCount = 0;
setfile = {};
numParsTot = 0; % initial value for total params
numEqsAuxsTot = 0; % initial value for total ics
ics_found  = ~ismember('IC',{parset(:).type});
pars_found = ~ismember('PAR',{parset(:).type});
tots_found = false;
allnames = {parset.name};

% get par and ic totals
while ~file_done & ~tots_found
    fline = fgetl(fid);
    setfile = {setfile{:}, fline};
    lineCount = lineCount + 1;
    if ~ischar(fline) % then EOF
        file_done = true;
        break
    end
    foundNumEqs = strfind(fline,'Number of equations and auxiliaries');
    if ~isempty(foundNumEqs)
        parsed = ParseParLine(fline);
        numEqsAuxsTot = parsed.num;
    end
    foundNumPars = strfind(fline,'Number of parameters');
    if ~isempty(foundNumPars)
        parsed = ParseParLine(fline);
        numParsTot = parsed.num;
    end
    tots_found = (numParsTot > 0 ) & (numEqsAuxsTot > 0);
end
if file_done % check this before totals
    disp(' Premature end of file reached. Cannot continue.')
    fclose(fid);
    return
end
if numEqsAuxsTot == 0 | numParsTot == 0 % then file is messed up
    disp(' The number of equations or parameters was not found in the SET file')
    fclose(fid);
    return
end

% ics come first in the SET file
while ~file_done & ~ics_found
    fline = fgetl(fid);
    setfile = {setfile{:}, fline};
    lineCount = lineCount + 1;
    if ~ischar(fline) % then EOF
        file_done = true;
        break
    end
    if strcmp(fline,'# Old ICs')
        ics_found = true;
    end
end
if file_done
    disp(' Premature end of file reached. ICs not found.')
    fclose(fid);
    return
end

% now change any specified ics, when encountered
icsDone = ~ics_found; % don't do the next code if no ics_found
icLineCount = 0;
if ~icsDone
	icnames = {};
	for parIx=1:lenpars % make a simple cell array list of parnames to be changed
        if strcmp(parset(parIx).type,'IC')
            icnames = {icnames{:},upper(parset(parIx).name)};
        end
	end
    numics = length(icnames);
	foundIcOccurrence = zeros(1,numics);
end

while ~icsDone & icLineCount <= numEqsAuxsTot
    fline = fgetl(fid);
    if ~ischar(fline) % then premature EOF
        disp(' Premature end of file reached while searching for initial conditions.')
        fclose(fid);
        return
    end
    icLineCount = icLineCount + 1;
    parsed = ParseParLine(fline);
    if parsed.num == NaN
        disp(' Problem parsing numerical value of parameter in SET file.')
        fclose(fid);
        return
    end
    [ispres ix] = ismember(parsed.name, upper(allnames));
    if ispres % then ic to be changed has been found
        if strcmp(parset(ix).type,'PAR')
            fprintf(' Found name `%s` that is in SET file as an initial condition\n',parset(ix).name)
            fclose(fid);
            return
        end
        foundIcOccurrence(ix) = 1;
        flineNew = [ num2str(parset(ix).val) '   ' parsed.name ];
    else
        flineNew = fline;
    end
    setfile = {setfile{:}, flineNew};
    lineCount = lineCount + 1;
    if sum(foundIcOccurrence) == numics % then made all changes necessary
        icsDone = true;
    end
end
if sum(foundIcOccurrence) ~= numics
    disp(' Error: Some i.c.`s specified for change were not found in SET file.')
    fclose(fid);
    return
end

%%%%%%

% now jump to parameter changes
while ~file_done & ~pars_found
    fline = fgetl(fid);
    if ~ischar(fline) % then EOF
        file_done = true;
        break
    end
    setfile = {setfile{:}, fline};
    lineCount = lineCount + 1;
    if strcmp(fline,'# Parameters')
        pars_found = true;
    end
end
if file_done
    disp(' Premature end of file reached. Parameters not found.')
    fclose(fid);
    return
end

% now change specified pars when encountered
parsDone = ~pars_found;
parLineCount = 0;
if ~parsDone
	parnames = {};
	for parIx=1:lenpars % make a simple cell array list of parnames to be changed
        if strcmp(parset(parIx).type,'PAR')
            parnames = {parnames{:},parset(parIx).name};
        end
	end
    numparams=length(parnames);
	foundParOccurrence = zeros(1,numparams);
end

while ~parsDone & parLineCount <= numParsTot
    fline = fgetl(fid);
    if ~ischar(fline) % then premature EOF
        disp(' Premature end of file reached while searching parameters.')
        fclose(fid);
        return
    end
    parLineCount = parLineCount + 1;
    parsed = ParseParLine(fline);
    if parsed.num == NaN
        disp(' Problem parsing numerical value of parameter in SET file.')
        fclose(fid);
        return
    end
    [ispres ix] = ismember(parsed.name, allnames);
    if ispres % then par to be changed has been found
        if strcmp(parset(ix).type,'IC')
            fprintf(' Found name `%s` that is in SET file as a parameter\n',parset(ix).name)
            fclose(fid);
            return
        end
        foundParOccurrence(ix) = 1;
        flineNew = [ num2str(parset(ix).val) '   ' parsed.name ];
    else
        flineNew = fline;
    end
    setfile = {setfile{:}, flineNew};
    lineCount = lineCount + 1;
    if sum(foundParOccurrence) == numparams % then made all changes necessary
        parsDone = true;
    end
end
if sum(foundParOccurrence) ~= numparams
    disp(' Error: Some parameters specified for change were not found in SET file.')
    fclose(fid);
    return
end

% file_done is still not true if we get to this point
while ~file_done
    fline = fgetl(fid);
    if ~ischar(fline) % then EOF
        file_done = true;
    else
        setfile = {setfile{:}, fline}; % finish reading in file
        lineCount = lineCount + 1;
    end
end
fclose(fid);

% write out changed file
fidw = fopen( filename ,'w');
if fidw ~= -1
	for line=1:lineCount
        fprintf(fidw, '%s\r\n', setfile{line});
	end
	fclose(fidw);
else
    fprintf(' Problem opening file %s for writing. Cannot continue\n',filename);
    return
end
success = true;
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function parsed = ParseParLine(fline)
parsed.num = NaN;
parsed.name = '';

[numStr nameStr] = strtok(fline);
if isempty(numStr) | isempty(nameStr)
    return
end
if ~isNum(numStr)
    return
end

nameStr = nameStr(isspace(nameStr)==0); % get rid of all whitespace

parsed.num = str2num(numStr);
parsed.name = nameStr;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
