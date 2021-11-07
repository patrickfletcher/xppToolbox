function xppdata=parseODEfile(filename)
%reads an xpp file, storing all information found. This version also
%accepts non-XPP range syntax introduced by Patrick Fletcher
%
% Output: xppdata, a struct containing the following structs:
%       par - parameter information (may be modified)
%           par.name, par.value, par.lb, par.ub
%           (XPP: p* name = val)
%       num - constant parameter information (hard coded as number)
%           num.name, num.value
%           (XPP: n* name=val)
%       var - variable information
%           var.name, var.value, var.formula, var.lb, var.ub
%           (XPP: name'=formula, dnamedt=formula; name(0)=value, init name=value)
%       fixed - fixed value (temporary quantity)
%           fixed.name, fixed.formula
%           (XPP: name=formula, !name = formula)
%       aux - user output
%           aux.name, aux.formula
%           (XPP: aux name=formula)
%
% A limited number of XPP options are also recognized. All options found
% (used in MatLab or not) are stored in: opt
%       opt - options
%           opt.name where name is an option name.
%           (XPP: @ name = value, value can be numeric or string)
%
% Comments found in the code are also collected.
%       comment - comments
%           comment.text
%           (XPP: # text, stores comments - combined successive comment lines)


% TODO
%
% Add help section about functions func(arg1,arg2,...)=formula
%
%  Add support for more XPP functionalities
%   - table - tables/interpolation [two versions: file, function]
%   - wiener - normally distributed random variable N(0,sqrt(dt))
%   - array/expansion of variables %[
%
%  CHECK all the reserved math and XPP words. Must either support them or
%  throw an error
%
%  Switch to using switch/case for different types of tokens (or enums)
%
%  Put things into better organization of subrountines to make main
%  function more clear
%
%

% Input check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if no filename was supplied, get one
if ~exist('filename','var')||isempty(filename)
    [name,path]=uigetfile('.ode','Select an ODE file');
    if ~ischar(name)
        error('File selection canceled, quitting...')
    end
    filename=fullfile(path,name);
end

% input file checking
if ~strcmp('.ode',filename(end-3:end))
    error(' Problem with filename: must end in .ode')
end
if exist(filename,'file') ~= 2
    error(' Problem with filename: file does not exist in specified path')
end

%try opening the file
fid = fopen(filename,'r');
if fid<0
    error([filename ' could not be opened.']);
end

[~,name]=fileparts(filename);

%storage structures
xppdata.name=name;

%these will each be fields in xppdata:
par=struct('name',{},'value',{},'lb',{},'ub',{});
num=struct('name',{},'value',{});
var=struct('name',{},'value',{},'lb',{},'ub',{},'formula',{},'formulaToken',{},'tokenType',{},'tokenIx',{});
fixed=struct('name',{},'formula',{},'formulaToken',{},'tokenType',{},'tokenIx',{});
func=struct('name',{},'full_arglist',{},'arg_names',{},'formula',{},'formulaToken',{},'tokenType',{},'tokenIx',{});
aux=struct('name',{},'formula',{},'formulaToken',{},'tokenType',{},'tokenIx',{});
wiener=struct('name',{});
opt=setXPPopt();
% opt=struct('name',{},'value',{});
comment=struct('text',{});

userNames={}; %set of all user defined names found
[ReservedNames,SimpleMathChars]=loadSpecialNames(); % Reserved keywords and math operations

file_done = false;
lineCount = 0; %current line
lastLineWasComment = false;

% Single pass, parse line by line. Semantic checking of formulas follows
while ~file_done
    
    fline = fgetl(fid);
    fullLine=fline;
    
    lineCount = lineCount + 1;
    
    
    %check to see if we reached the end of the ODE file yet
    if ~ischar(fline)
        %         disp(['Finished parsing ' filename ' - now checking formulas... '])
        file_done=true;
        break
    end
    
    %empty lines are ignored, non-empty lines have leading/trailing spaces removed
    fline = strtrim(fline);
    if isempty(fline)
        lastLineWasComment=false;
        %skip to next line if blank
        continue
    end
    
    %%% Lines with a unique first character, no space required after.
    %
    % # - comment
    % " - comment (displayable in XPP, potentially with actions)
    % !name=formula - derived parameter
    % 0=expr - algebraic equation
    
    if strcmp(fline(1),{'#'}) %comment
        if lastLineWasComment
            comment(end).text=[comment(end).text ' ' fline(2:end)];
        else
            comment(end+1).text=fline(2:end);
        end
        
        lastLineWasComment=true;
        %nothing more of interest in this line
        continue
        
    elseif strcmp(fline(1),{'"'}) %comment
        comment(end+1).text=fline(2:end);
        
        lastLineWasComment=true;
        %nothing more of interest in this line
        continue
        
    elseif strcmp(fline(1),{'%'}) % && ~strcmp(fline(2),{'['}) %comment
        comment(end+1).text=fline(2:end);
        
        lastLineWasComment=true;
        %nothing more of interest in this line
        continue
        
    elseif strcmp(fline(1),{'$'}) %comment??? seen this a few times...
        comment(end+1).text=fline(2:end);
        
        lastLineWasComment=true;
        %nothing more of interest in this line
        continue
        
    elseif strcmpi(fline(1), '!')
        [LHS, RHS] = ExtractFormula(lower(fline(2:end)));
        %assume LHS is a name
        if isValidName(LHS)
            fixed(end+1).name=LHS;  % Treat as a fixed quantity for now
            fixed(end).formula=RHS;
            userNames{end+1}=LHS;
        end
        lastLineWasComment=false;
        continue
        
    elseif strcmpi(fline(1), '0') %not yet supported
        disp([num2str(lineCount) ': ' fullLine])
        error('Algebraic equations not yet supported')
    end
    
    %getting this far means this line isn't blank or a comment
    lastLineWasComment=false;
    
    %XPP is case insensitive, force all lower.
    fline = lower(fline);
    
    %multiline expansion with %[j1..j2] syntax
    if length(fline)>1 && strcmp(fline(1:2),'%[')
        disp([num2str(lineCount) ': ' fullLine])
        error('Variable expansion with %[j1..j2] ... % is not yet supported')
        
        %this option would require lookahead in lines to find the closing %
    end
    
    % KEYWORD lines (have a special word/char as the first word. Space
    % required between keyword and following info.
    %
    % @ name=value - options
    % a name=formula - auxiliary variable
    % b - boundary condition
    % done - end of file  
    % e - export
    % g - global
    % i name=value - init 
    % m - markov
    % n name=value - number (constant parameter)
    % o - option filename (deprecated)
    % p name=value - parameter
    % se - set (specify a named set of param/init/option)
    % so - solve algebraic condition
    % sp - special function
    % t - table
    % v - volterra
    % w name - wiener
    
    [token, rest] = strtok(fline); %token is the word left of the space
    
    token=strtrim(token); %remove leading/trailing spaces (Can I remove all spaces?)
    rest=strtrim(rest);
    
    %Make sure it is a keyword line
    if ~isempty(rest) && ~ismember('=', token) && ~strcmpi(rest(1),'=')
        
        % @ - options
        if strcmpi(token(1),'@')
            [parsed, numparsed] = ParseLine(rest,2,lineCount);
            for i=1:numparsed
                opt=setXPPopt(opt,parsed(i).name, parsed(i).num);
%                 opt(end+1).name=parsed(i).name;
%                 opt(end).value=parsed(i).num;
            end
            
            continue
            
            % a - auxiliary variable
        elseif strcmpi(token(1), 'a')
            [LHS,RHS]=ExtractFormula(rest);
            aux(end+1).name=LHS;
            aux(end).formula=RHS;
            
            continue
            
            % b - boundary condition
        elseif strcmpi(token(1), 'b')
            disp([num2str(lineCount) ': ' fullLine])
            error('Boundary Condition is not yet supported')
            
        elseif strcmpi(token, 'done')
            
            file_done=true;
            break;
            
            % e - export
        elseif strcmpi(token(1), 'e')
            disp([num2str(lineCount) ': ' fullLine])
            error('Export is not yet supported')
            
            % g - global
        elseif strcmpi(token(1), 'g')
            disp([num2str(lineCount) ': ' fullLine])
            error('Global is not yet supported')
            
            % i - init
        elseif strcmpi(token(1), 'i')
            [parsed, numparsed] = ParseLine(rest,1,lineCount);
            for i=1:numparsed
                %Check if this var is found. If not, add it to the list
                [varFound,varLoc]=ismember(parsed(i).name,{var(:).name});
                if varFound
                    var(varLoc).value=parsed(i).num;
                    var(varLoc).lb=parsed(i).lb; %empty if no range specified
                    var(varLoc).ub=parsed(i).ub;
                else
                    if isValidName(parsed(i).name)
                        var(end+1).name=parsed(i).name;
                        userNames{end+1}=parsed(i).name;
                        var(end).value=parsed(i).num;
                        var(end).lb=parsed(i).lb; %empty if no range specified
                        var(end).ub=parsed(i).ub;
                    end
                end
            end
            
            continue
            
            % m - markov
        elseif strcmpi(token(1), 'm')
            disp([num2str(lineCount) ': ' fullLine])
            error('Markov is not yet supported')
            
            % n - number (constant parameter)
        elseif strcmpi(token(1), 'n')
            [parsed, numparsed] = ParseLine(rest,0,lineCount);
            for i=1:numparsed
                if isValidName(parsed(i).name)
                    num(end+1).name=parsed(i).name;
                    userNames{end+1}=parsed(i).name;
                    num(end).value=parsed(i).num;
                end
            end
            
            continue
            
            % o - option filename (deprecated)
        elseif strcmpi(token(1), 'o')
            disp([num2str(lineCount) ': ' fullLine])
            error('Option Filename is not supported')
            
            % p - parameter
        elseif strcmpi(token(1), 'p')
            [parsed, numparsed] = ParseLine(rest,0,lineCount);
            for i=1:numparsed
                if isValidName(parsed(i).name)
                    par(end+1).name=parsed(i).name;
                    userNames{end+1}=parsed(i).name;
                    par(end).value=parsed(i).num;
                    par(end).lb=parsed(i).lb;
                    par(end).ub=parsed(i).ub;
                    
                end
            end
            
            continue
            
            % se - set (specify a named set of param/init/option)
        elseif length(token)>1 && strcmpi(token(1:2), 'se')
            disp([num2str(lineCount) ': ' fullLine])
            error('Set is not yet supported')
            
            % so - solve algebraic condition
        elseif length(token)>1 && strcmpi(token(1:2), 'so')
            disp([num2str(lineCount) ': ' fullLine])
            error('Solve is not yet supported')
            
            % sp - special function
        elseif length(token)>1 && strcmpi(token(1:2), 'sp')
            disp([num2str(lineCount) ': ' fullLine])
            error('Special functions are not yet supported')
            
            % t - table
        elseif strcmpi(token(1), 't')
            disp([num2str(lineCount) ': ' fullLine])
            error('Tables are not yet supported')
            
            % v - volterra
        elseif strcmpi(token(1), 'v')
            disp([num2str(lineCount) ': ' fullLine])
            error('Volterra equations are not yet supported')
            
            % w - wiener
        elseif strcmpi(token(1), 'w')
            [parsed, numparsed] = ParseLine(rest,3,lineCount);
            for i=1:numparsed
                if isValidName(parsed(i).name)
                    wiener(end+1).name=parsed(i).name;
                    userNames{end+1}=parsed(i).name;    
                end
            end
            
            continue
            
        else
            disp([num2str(lineCount) ': ' fullLine])
            error(['The keyword "' token '" is unrecognized. Cannot continue'])
        end
    end
    
    if strcmpi(token, 'done')
        file_done=true;
        break;
    end
    
    %if no keywords were found, the only remaining possibility is a variant
    %of 'name=formula'. name can represent either a variable or a fixed
    %quantity.
    
    %%% NAME=FORMULA variants:
    % dname/dt=formula            - Differential equation (defines var)
    % name'=formula               - Differential equation (defines var)
    % name(0)=value/expr          - initial data. value defines var, expr=>DDE init.
    % name(t+1)=formula           - Difference equation (defines var)
    % name(t)=formula             - Volterra integral equation
    % name(var1,var2,...)=formula - function definition (var1-var9 names don't matter)
    % name[j1..j2]=expr[j]        - expansion with indexer [j]
    % name=formula                - temporary quantity (defines fixed)
    
    
    [LHS,RHS]=ExtractFormula(fline);
    
    if isempty(LHS) || isempty(RHS)
        disp([num2str(lineCount) ': ' fullLine])
        error('There was an empty left or right hand side')
    end
    
    isInit=false;
    isEqn=false;
    isFixed=false;
    isFunc=false; arg_names=[];
    
    %check if LHS indicates ODE, IC, or fixed value
    
    %dname/dt=formula
    if length(LHS)>4 && strcmpi(LHS(1),'d') && strcmpi(LHS(end-2:end),'/dt')
        isEqn=true;
        name=LHS(2:end-3);
        
        %name'=formula
    elseif length(LHS)>1 && LHS(end)==''''
        isEqn=true;
        name=LHS(1:end-1);
        
        %name(0)=formula
    elseif length(LHS)>3 && strcmpi(LHS(end-2:end),'(0)')
        isInit=true;
        name=LHS(1:end-3);
        
        %name(t+1)=formula
    elseif length(LHS)>5 && strcmpi(LHS(end-4:end),'(t+1)')
        disp([num2str(lineCount) ': ' fullLine])
        error('Difference equations are not yet supported');
        
        %name(t)=formula
    elseif  length(LHS)>3 && strcmpi(LHS(end-2:end),'(t)')
        disp([num2str(lineCount) ': ' fullLine])
        error('Volterra integral equations are not yet supported');
        
        %name(var1,var2,...)=formula
    elseif length(LHS)>3 && any(LHS=='(') && LHS(end)==')'
        isFunc=true;
        ixLeftPar=find(LHS=='(');
        name=LHS(1:ixLeftPar-1);
        
        full_arglist=LHS(ixLeftPar+1:end-1); %includes the commas, useful for forming MatLab anonymous function.
        
        rest=LHS(ixLeftPar+1:end-1); %now grab the names of each argument.
        while ~isempty(rest)
            [arg_names{end+1},rest]=strtok(rest,',');
        end
        
        %name[j1..j2]=expr[j]
    elseif  length(LHS)>5 && any(LHS=='[')
        disp([num2str(lineCount) ': ' fullLine])
        error('Expansions using name[j1..j2]=expr[j] are not yet supported');
        
        %name=formula
    else
        isFixed=true;
        name=LHS;
    end
    
    %%% now process the RHS
    
    if isInit
        
        if any(RHS=='[')
            [value,lb,ub,restStr]=parseRangeToken(RHS);
        elseif isNumericLiteral(RHS)
            value=str2double(RHS);
            restStr=[];
            lb=[];
            ub=[];
        else
            disp([num2str(lineCount) ': ' fullLine])
            error('I don''t understand what you put in the formula above!')
        end
        
        if ~isempty(restStr)
            disp([num2str(lineCount) ': ' fullLine])
            error('Only one initial condition with syntax name(0)=value allowed per line')
        end
        
        % check if we already found the variable
        [varFound,varLoc]=ismember(name,{var(:).name});
        if varFound
            var(varLoc).value=value;
            var(varLoc).lb=lb;
            var(varLoc).ub=ub;
        else
            if isValidName(name)
                %if not a name used already as a param/fixed/reserved
                var(end+1).name=name;
                userNames{end+1}=name;
                var(end).value=value;
                var(end).lb=lb;
                var(end).ub=ub;
            end
        end
        
    elseif isEqn
        % check if we already found the variable
        [varFound,varLoc]=ismember(name,{var(:).name});
        if varFound
            var(varLoc).formula=RHS;
        else
            if isValidName(name)
                var(end+1).name=name;
                userNames{end+1}=name;
                var(end).formula=RHS;
            end
        end
        
    elseif isFixed
        if isValidName(name)
            fixed(end+1).name=name;
            userNames{end+1}=name;
            fixed(end).formula=RHS;
        end
        
    elseif isFunc
        if isValidName(name)
            func(end+1).name=name;
            userNames{end+1}=name;
            func(end).arg_names=arg_names; %for checking validity of formula
            func(end).full_arglist=full_arglist;
            func(end).formula=RHS;
        end
    end
end

fclose(fid);

% Extract location and type information of each of the pieces of the
% formula and store them (for conversions to other languages). At the same
% time make sure each formula uses only valid names.

% (tokenization of formulae and semantic checking)

% TODO: proper check that formulas make mathematical sense!

nPar=length(par);
nNum=length(num);
nVar=length(var);
nFixed=length(fixed);
nFunc=length(func);
nAux=length(aux);
nWiener=length(wiener);
nOpt=length(opt);

for i=1:nPar
    if isempty(par(i).lb)
        par(i).lb=par(i).value;
    end
    if isempty(par(i).ub)
        par(i).ub=par(i).value;
    end
end


i=1;
while i<=nFixed
    thisName=fixed(i).name;
    thisFormula=fixed(i).formula;
    [tokens,types,ix]=TokenizeFromula(thisName,thisFormula);
    
    %reorder fixed formulas based on dependencies of the right hand side
    fixDeps=ix(types==4);
    if ~isempty(fixDeps)
        maxDep=max(fixDeps);
        if maxDep>i
            fixed=[fixed(1:i-1),fixed(maxDep),fixed(i:maxDep-1), fixed(maxDep+1:end)]; %move the worst dependency up 
%             fixed=[fixed(1:i-1), fixed(i+1:maxDep), fixed(i), fixed(maxDep+1:end)]; %move this below its latest depenency
            continue
        end
    end
    
    fixed(i).formulaToken=tokens;
    fixed(i).tokenType=types;
    fixed(i).tokenIx=ix;
    
    i=i+1;
end


for i=1:nVar
    thisName=var(i).name;
    thisFormula=var(i).formula;
    [tokens,types,ix]=TokenizeFromula(thisName,thisFormula);
    var(i).formulaToken=tokens;
    var(i).tokenType=types;
    var(i).tokenIx=ix;
    if isempty(var(i).value) %if no Initial value, default to 0 like xpp does
        var(i).value=0;
    end
    if isempty(var(i).lb)
        var(i).lb=var(i).value;
    end
    if isempty(var(i).ub)
        var(i).ub=var(i).value;
    end
end


for i=1:nFunc
    % XPP doesn't care what the args are called, as long as rhs tokens are
    % found names or elements of arglist. Here we will rename the arguments
    % to formal names arg1-arg9 (like XPP) to avoid accidental naming
    % problems (e.g. argname=y or p would crash .m file, since vector of
    % state variables and parameters are called y and p)
    
    thisName=func(i).name;
    thisFormula=func(i).formula;
    [tokens,types,ix]=TokenizeFromula(thisName,thisFormula,func(i).arg_names);
    func(i).formulaToken=tokens;
    func(i).tokenType=types;
    func(i).tokenIx=ix;
end

for i=1:nAux
    thisName=aux(i).name;
    thisFormula=aux(i).formula;
    [tokens,types,ix]=TokenizeFromula(thisName,thisFormula);
    aux(i).formulaToken=tokens;
    aux(i).tokenType=types;
    aux(i).tokenIx=ix;
end


%all checks successful: package the output
xppdata.par=par;
xppdata.num=num;
xppdata.var=var;
xppdata.fixed=fixed;
xppdata.func=func;
xppdata.aux=aux;
xppdata.opt=opt;
xppdata.comment=comment;

%convenience outputs
xppdata.nPar=nPar;
xppdata.parNames={par(:).name};
xppdata.p0=[par(:).value];

xppdata.nNum=nNum;
xppdata.numNames={num(:).name};
xppdata.n0=[num(:).value];

xppdata.nVar=nVar;
xppdata.varNames={var(:).name};
xppdata.x0=[var(:).value];

xppdata.nFixed=nFixed;
xppdata.fixNames={fixed(:).name};

xppdata.nFunc=nFunc;
xppdata.funcNames={func(:).name};

xppdata.nAux=nAux;
xppdata.auxNames={aux(:).name};

xppdata.nWiener=nWiener;
xppdata.WienerNames={wiener(:).name};

xppdata.nOpt=nOpt;
xppdata.optNames=fieldnames(opt);
% xppdata.optVals=[opt(:).value]; %

%if we got to the end, then all formulas are valid!
disp('All formulas are valid! Parsing Successful.')


disp('')

% success=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function result=isValidName(name)
        if any(strcmpi(name,userNames))
            disp([num2str(lineCount) ': ' fullLine])
            error('Error parsing ODE file: duplicate name')
        elseif any(strcmpi(name,ReservedNames))
            disp([num2str(lineCount) ': ' fullLine])
            error('Error parsing ODE file: illegal name')
        end
        result=true;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [tokens,types,ix] = TokenizeFromula(name,formula,arglist)
        %checks formula for internal consistency and then returns a list of
        %tokens, their types, and the index into list of all found elements
        %of that type.
        
        %types:
        %0=number, 1=simplemath, 2=mathname, 3=par, 4=fixed, 5=var,
        %6=function, 7=function argument, 8=wiener
        
        %ix is index into the array of given type
        
        if ~exist('arglist','var'), arglist={};end
        
        full_formula=formula; %for error message, if needed
        formula = formula(isspace(formula)==0); % get rid of whitespace
        
        %regex codes:
        real_in = '\d+';
        real_dn = '(\d+\.\d*|\d*\.\d+)';
        real_sn = '(\d*\.?\d+[Ee][+\-]?\d+|\d+\.?\d*[Ee][+\-]?\d+)';
        real_num = ['(', real_sn, '|', real_dn, '|', real_in, ')']; %-? I let the unary "-" fall through to math ops...
        word = '[a-zA-Z]+\w*';
        ops = '[()^*/+-<>]';
        
        %tokenize
        tokens = regexp(formula,[real_num,'|',word,'|',ops],'match');
        
        %parentheses check
        numLeftParenth=nnz(formula=='(');
        numRightParenth=nnz(formula==')');
        if numLeftParenth~=numRightParenth
            disp([name,formula])
            error('Parentheses don''t match in the formula above');
        end
        
        types=zeros(size(tokens));
        ix=zeros(size(tokens));
        
        %read the line and parse out all tokens
        for j=1:length(tokens)
            
            tok = tokens{j};
            %check the token type and if applicable find its index
            if isNumericLiteral(tok) 
                tokens{j}=str2double(tokens{j}); %convert string to a double
                types(j)=0;
                
            %is the token a math operation?
            elseif ismember(tok,SimpleMathChars) %one of (,),-,+,*,/,^,<,>
                types(j)=1;
                
            elseif regexp(tok,word)
                %is the word a reserved name? (includes math functions)
                [isTrue,thisIx]=ismember(tok,ReservedNames);
                if isTrue
                    types(j)=2;
                    ix(j)=thisIx;
                    continue; %skip to next token. Alt: check all to ensure no re-use?
                end
                
                %is the word a named constant?
                [isTrue,thisIx]=ismember(tok,{num(:).name});
                if isTrue
                    types(j)=0;
                    ix(j)=thisIx;
                    continue;
                end
                   
                %is the word a named parameter
                [isTrue,thisIx]=ismember(tok,{par(:).name});
                if isTrue
                    types(j)=3;
                    ix(j)=thisIx;
                    continue;
                end
                      
                %is the word a named formula?
                [isTrue,thisIx]=ismember(tok,{fixed(:).name});
                if isTrue
                    types(j)=4;
                    ix(j)=thisIx;
                    continue;
                end

                %is the word a named variable
                [isTrue,thisIx]=ismember(tok,{var(:).name});
                if isTrue
                    types(j)=5;
                    ix(j)=thisIx;
                    continue;
                end

                %is the word a named function?
                [isTrue,thisIx]=ismember(tok,{func(:).name});
                if isTrue
                    types(j)=6;
                    ix(j)=thisIx;
                    continue;
                end

                %is the word a named function argument?
                [isTrue,thisIx]=ismember(tok,arglist);
                if isTrue
                    types(j)=7;
                    ix(j)=thisIx;
                    continue;
                end

                %is the word a named Wiener variable?
                [isTrue,thisIx]=ismember(tok,{wiener(:).name});
                if isTrue
                    types(j)=8;
                    ix(j)=thisIx;
                    continue;
                end
                
                %no other cases are defined yet
            else
                %if we get this far, the name is not one we've found yet.
                disp(['---> ' tok ' in line ' name '=' full_formula ]);
                error('Error parsing ODE file: unknown name')
                
            end
        end
        
    end



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more helper functions, not nested
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [LHS, RHS] = ExtractFormula(string)
    idx_Eq=find(string=='=');

    if length(idx_Eq)>1
        disp(string)
        error(['Too many "=" in ''' string '''' ])
    end

    LHS=string(1:idx_Eq-1);
    RHS=string(idx_Eq+1:end);

    LHS=LHS(~isspace(LHS));  %now we can safely remove spaces
    RHS=RHS(~isspace(RHS));
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [parsed, numparsed] = ParseLine(full_line, type, lineCount)
    %type 0=par, 1=init, 2=opt, 3=wiener

    parseddef.num = NaN;
    parseddef.name = '';
    parseddef.lb = [];
    parseddef.ub = [];
    parsed=parseddef;
    thisL = full_line(isspace(full_line)==0); % get rid of whitespace
    thisL = strrep(thisL,'=',' '); % convert = to space
    thisL = strrep(thisL,',',' '); % convert , to space
    numparsed = 0;
    endofline = false;

    thisRest=thisL;

    while ~endofline
        [nameStr, thisRest] = strtok(thisRest);
        if isempty(nameStr)
            endofline = true;
            continue
        end

        numparsed = numparsed + 1;

        nextTok = strtok(thisRest);
        rangeBrackets=any(nextTok=='[');
        if rangeBrackets
            [numericValue,thislb,thisub,thisRest]=parseRangeToken(thisRest);
        else
            [numStr, thisRest] = strtok(thisRest);
            if type~=3 && isempty(numStr) || isempty(nameStr)
                disp([num2str(lineCount) ': ' full_line]);
                error('Param parse error: Parameter name or value missing!')
            end
            if type==0 && ~isNumericLiteral(numStr(isspace(numStr)==0))
                disp([num2str(lineCount) ': ' full_line]);
                error('Param parse error: Parameter value not a number!')
            end
            if type==1 && ~isNumericLiteral(numStr(isspace(numStr)==0))
                disp([num2str(lineCount) ': ' full_line]);
                error('Init parse error: Initial value not a number!')
            end
            if  type<2
                numericValue=str2double(numStr);
            else
                numericValue=numStr; %options are allowed to have char values
            end
        end

        parsed(numparsed).name = nameStr(isspace(nameStr)==0); % get rid of any whitespace
        parsed(numparsed).num  = numericValue;
        if rangeBrackets
            parsed(numparsed).lb = thislb;
            parsed(numparsed).ub = thisub;
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [num,lb,ub,restStr]=parseRangeToken(tok)
%this helper extracts # and range from token of form #[#,#]

tok = strtrim(tok); % get rid of whitespace
ixL=find(tok=='[',1,'first');
ixR=find(tok==']',1,'first');
numStr=tok(1:ixL-1);
if isNumericLiteral(numStr)
    num=str2double(numStr);
else
    disp(tok);
    error('Incorrect format for value range: name=value[lo,hi]')
end
range  = eval(tok(ixL:ixR));  % currently no error checking on this...
lb=range(1);
ub=range(2);
restStr=tok(ixR+1:end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = isNumericLiteral(str)
% parameter `data` is a string
str = string(str);

real_in = '\d+';
real_dn = '(\d+\.\d*|\d*\.\d+)';
real_sn = '(\d*\.?\d+[Ee][+\-]?\d+|\d+\.?\d*[Ee][+\-]?\d+)';
real_num = ['-?(', real_sn, '|', real_dn, '|', real_in, ')']; 

match = regexp(str,real_num,'match');

result = false;
if match == str
    result = true;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ReservedNames,SimpleMathChars,ReservedSymbols,FunctionsOneArg,FunctionsTwoArg,SpecialXPPSyntax]=loadSpecialNames()
%XPP keywords and built in math names that can't be used as user names.

%could create array of structs like Bard's, particularly for function
%names and operations, so we can check for the correct usage. (e.g. correct
%number of function arguments; if matches then,else; etc.)

%list of math and XPP special names
ReservedXPP = {'sin', 'cos', 'tan',...
    'asin', 'acos', 'atan', 'atan2', ...
    'sinh', 'cosh', 'tanh',...
    'exp', 'ln', 'log', 'log10',...
    'sqrt', 'heav', 'sign',...
    'max', 'min', 'mod', 'abs', 'ceil', 'flr',...
    't', 'pi', ...
    'delay', 'del_shft', 'shift', 'ishift',...
    'if', 'then', 'else', ...
    'ran', 'normal','lgamma', ...
    'besseli', 'besselj', 'bessely', 'erf', 'erfc', ...
    'arg1', 'arg2', 'arg3', 'arg4', 'arg5','arg6', 'arg7', 'arg8', 'arg9','arg10', ...
    'arg11', 'arg12', 'arg13', 'arg14','arg15', 'arg16', 'arg17', 'arg18','arg19','arg20',...
    '@', '$', '**', '~', ...
    '[',']','{','}',...
    '\#', 'int',...
    'sum', 'of','i'''...
    'hom_bcs',...
    };

%The following partition may come in handy for error checking

ReservedSymbols={'t', 'pi', ...
    '@', '$', '**', '~', ...
    'arg1', 'arg2', 'arg3', 'arg4', 'arg5','arg6', 'arg7', 'arg8', 'arg9','arg10', ...
    'arg11', 'arg12', 'arg13', 'arg14','arg15', 'arg16', 'arg17', 'arg18','arg19','arg20',...
    };


FunctionsOneArg={'sin', 'cos', 'tan',...
    'asin', 'acos', 'atan',...
    'sinh', 'cosh', 'tanh',...
    'exp', 'ln', 'log', 'log10',...
    'sqrt', 'heav', 'sign',...
    'abs', 'ceil', 'flr',...
    'ran',...
    'lgamma', 'erf', 'erfc', ...
    'hom_bcs',...
    };

FunctionsTwoArg={'atan2', ...
    'max', 'min', 'mod', ...
    'delay',...
    'normal',...
    'besseli', 'besselj', 'bessely',...
    'shift', 'ishift',...
    };

SpecialXPPSyntax={'del_shft',...
    'if', 'then', 'else', ...
    'sum', 'of','i'''...
    'int','[',']', '{','\#','}',...
    };


RelationalOperations={'|', '>', '<', '==', '>=', '<=', '!=', 'not','&'}; %what about double &&, ||

SimpleMathChars = {'(',')','^','*','/','+','-',','}; %comma for function arg separator

ReservedNames=[ReservedXPP, RelationalOperations];
end


function XPPopt=setXPPopt(XPPopt,name,value)
% Need to figure out how sloppy XPP is with naming options, so I can
% emulate that here. Perhaps only the first few letters are needed...

if ~exist('XPPopt','var') || isempty(XPPopt)
    XPPopt=struct('allopt',struct('name',{},'value',{}),...
        'method',[],...
        'dt',0.05, 'dtmin',1e-6, 'dtmax',10,'toler',1e-3,'atoler',1e-6,...
        't0',0, 'trans',0, 'total',150,...
        'maxstor',100000,'bound',100,'nout',1,...
        'seed',[],...
        'nplot',1,...
        'xp','t','yp',[],'zp',[],...'nplot',1,...'xp2',[],'yp2',[],'zp2',[],... %do more of these? how does XPP do it?
        'xlo',[],'xhi',[],'ylo',[],'yhi',[],...
        'axes',2,'phi',[],'theta',[],...
        'xmin',[],'xmax',[],'ymin',[],'ymax',[],'zmin',[],'zmax',[],...
        'fractionYup',0.35,'fractionYdown',0.25,'fractionDYup',0.0,'fractionDYdown',0.0,...%<--mine, not xpp options!
        'minYamp',0.0,'minDYamp',0.0,'minIMI',0.0,...%<--mine, not xpp options!
        'normTol',0.1); %<--mine, not xpp options!
end

if exist('name','var')&&~isempty(name)&&exist('value','var')&&~isempty(value)
    switch lower(name)
        case {'meth','method'}
            XPPopt.method=value;
        case 'dt'
            XPPopt.dt=str2double(value);
        case 'dtmin'
            XPPopt.dtmin=str2double(value);
        case 'dtmax'
            XPPopt.dtmax=str2double(value);
        case 'toler'
            XPPopt.toler=str2double(value);
        case 'atoler'
            XPPopt.atoler=str2double(value);
        case 't0'
            XPPopt.t0=str2double(value);
        case 'trans'
            XPPopt.trans=str2double(value);
        case 'total'
            XPPopt.total=str2double(value);
        case 'maxstor'
            XPPopt.maxstor=str2double(value);
        case {'bound','bounds'}
            XPPopt.bound=str2double(value);
        case {'nout','njmp'}
            XPPopt.nout=str2double(value);
        case 'seed'
            XPPopt.seed=str2double(value);
        case 'nplot'
            XPPopt.nplot=str2double(value);
        case 'xp'
            XPPopt.xp=value;
        case 'yp'
            XPPopt.yp=value;
        case 'zp'
            XPPopt.zp=value;
        case 'xlo'
            XPPopt.xlo=str2double(value);
        case 'xhi'
            XPPopt.xhi=str2double(value);
        case 'ylo'
            XPPopt.ylo=str2double(value);
        case 'yhi'
            XPPopt.yhi=str2double(value);
        case 'axes'
            XPPopt.axes=str2double(value);
        case 'phi'
            XPPopt.phi=str2double(value);
        case 'theta'
            XPPopt.theta=str2double(value);
        case 'xmin'
            XPPopt.xmin=str2double(value);
        case 'xmax'
            XPPopt.xmax=str2double(value);
        case 'ymmin'
            XPPopt.ymin=str2double(value);
        case 'ymax'
            XPPopt.ymax=str2double(value);
        case 'zmin'
            XPPopt.zmin=str2double(value);
        case 'zmax'
            XPPopt.zmax=str2double(value);
            
            %now my options
        case 'fractionyup'
            XPPopt.fractionYup=str2double(value);
        case 'fractionydown'
            XPPopt.fractionYdown=str2double(value);
        case 'fractiondyup'
            XPPopt.fractionDYup=str2double(value);
        case 'fractiondydown'
            XPPopt.fractionDYdown=str2double(value);
        case 'minyamp'
            XPPopt.minYamp=str2double(value);
        case 'mindyamp'
            XPPopt.minDYamp=str2double(value);
        case 'minimi'
            XPPopt.minIMI=str2double(value);
        case 'normtol'
            XPPopt.normTol=str2double(value);
            
        otherwise
            XPPopt.allopt(end+1).name=name;
            XPPopt.allopt(end).value=value;
    end
end
end