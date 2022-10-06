function [destinationFile, xppdata]=ode2m(source, destinationFile,options)
arguments
    source=[]
    destinationFile=[]
    options.verbose=true
    options.vectorized=true
end
% code conversion from ODE file for XPP to a .m file that can be run with
% MatLab. Calls parseODEfile, writes the ODE information to a .m file,
% then packages the secondary parameters into a vector pars. Also produced
% is a structure containing options found in the ODE file (e.g. dt, total,
% clODE stuff, etc.)
%
% USAGE
%   [mFunctionName, xppdata]=ode2m(arg1, destinationFile)
%
% Inputs
%   none - will prompt for ODE file, write to same name .m in working dir
%   source: either full path to ODE file, or pre-generated xppdata struct
%           from parseODEfile 
%   destinationFile: specified output filename. if omitted, a default will be generated based on the
%                    source filename
%   verbose {true}/false - toggle display of information in the console
%   vectorized {true}/false - generates a function compatible with
%                             "vectorized" option in Matlab ODE solvers
%                             (internally, they use rows for variables)
%
% Outputs
%   mFileName: string conatining the destination filename without extension
%   xppdata: structure containing all ODE file information
%
% Example for using MatLab's ODE solvers:
%   [mFunctionName,p,y0]=ode2m('./xppSrc/lactotroph.ode');
%   fun=eval(['@(t,y)' mFunctionName '(t,y,p)']);
%   ode45(fun,[0,1000],y0);
%
%



%TO DO
% - destination file overwriting behavior? now destroys any existing copy
% - don't use symbolic toolbox heaviside!!

verbose=options.verbose;

fileExtension='m';

% input checking
if nargin==0 || isempty(source)
    [name,path]=uigetfile('.ode','Select an ODE file');
    if ~ischar(name)
        disp('File selection canceled, quitting...')
        return
    end
    source=fullfile(path,name);
end

if ischar(source) %filename input: parse it.
    srcfilename=source;
    % Extract info from the ODE file
    xppdata=parseODEfile(srcfilename);
else %assume it is ivpStruct (TODO: errorchecking)
    xppdata=source;
end

xppdata.fun=eval("@(t,x)"+xppdata.name+"(t,x,xppdata.p0)");

%Build the destination filename
if isempty(destinationFile)
    destinationFile=xppdata.name;
else
    %destinationFile was provided. Extract just the filename without
    %extension or path.
    [~,destinationFile]=fileparts(destinationFile);
    
    %check whether it exists, and ask if overwriting it is ok?
end

output_file=BuildOutputFile(xppdata,destinationFile,options);
lineCount=length(output_file);

%delete a previous version - silently destroys any previous version!
fullPath=[pwd filesep destinationFile '.' fileExtension];
if exist(fullPath,'file')==2
    delete(fullPath)
end

% write to file
fidw=fopen(fullPath,'w');

if fidw > -1
    for line=1:lineCount
        fprintf(fidw, '%s\n', output_file{line});
    end
    status=fclose(fidw);
else
    fprintf(' Problem opening file %s for writing. Cannot continue\n',fullPath);
    return
end


%check to make sure file is written
ready=false;
while ~ready
    status=exist(fullPath,'file');
    if status
        ready=true;
    end
end

%summary information
if verbose
    disp('Found:')
    disp([num2str(xppdata.nVar) ' variables'])
    disp([num2str(xppdata.nPar) ' parameters'])
    disp([num2str(xppdata.nNum) ' numbers'])
    disp([num2str(xppdata.nFixed) ' fixed quantities'])
    disp([num2str(xppdata.nFunc) ' functions'])
    disp([num2str(xppdata.nWiener) ' Wiener variables'])
    disp([num2str(xppdata.nAux) ' auxiliary variables'])
    
    %give suggestions on which matlab solver to try based on the method option
    %detected. TODO: expand this with a better understanding of exactly which
    %solvers XPP is using!
    methodHint='';
    xppmethod=xppdata.opt.method;
    if xppdata.nWiener>0
        methodHint=['Wiener variables detected. Use ode_euler'];
%     elseif ~isempty(xppdata.opt.method)
%         methodHint=['XPP method option was set to ' xppdata.opt.method '. '];
%         switch xppdata.opt.method
    elseif ~isempty(xppmethod)
        methodHint=['XPP method option was set to ' xppmethod '. '];
        switch xppmethod
            case {'euler'}
                methodHint=[methodHint, 'Try ode_euler'];
            case {'runge','rungekutta','modeuler','5dp','83dp'}
                methodHint=[methodHint, 'Try ode45'];
            case {'backeuler','gear','stiff','cvode'}
                methodHint=[methodHint, 'Try ode15s'];
            case '2rb'
                methodHint=[methodHint, 'Try ode23s'];
        end
    end
    disp(methodHint)
    
end

end

function output_file=BuildOutputFile(xppdata,destinationFile,options)

num=xppdata.num;
var=xppdata.var;
fixed=xppdata.fixed;
func=xppdata.func;
aux=xppdata.aux;

% build the output as a cell array containing each line.
% Structure:
%   1. function signature
%   2. initialize output variables
%   3. anonymous function declaration
%   4. fixed formulas
%   5. differential equations
%   6. auxiliary output variables

%names for MatLab internal vectors, to try to avoid potential name collisions
%with user fixed quantities, parameters, or variables.
stateName='x_';
paramName='p_';
slopeName='dx_';
auxName='aux_';
wienerName='w_';


if xppdata.nWiener==0
output_file{1}=['function [' slopeName ', ' auxName ']=' destinationFile '(t,' stateName ',' paramName ')'];
else
output_file{1}=['function [' slopeName ', ' auxName ']=' destinationFile '(t,' stateName ',' wienerName ',' paramName ')'];
end

%if vectorized, vars need to be rows
output_file{end+1}='if any(size(x_)==1), x_=x_(:); end';
output_file{end+1}=[slopeName '=zeros(' num2str(xppdata.nVar) ',size(x_,2));'];
output_file{end+1}=[auxName '=zeros(' num2str(xppdata.nAux) ',size(x_,2));'];


%put anonymous function definitions first
for i=1:xppdata.nFunc
    
    name=func(i).name; %case is all lower after parser (XPP is case insensitive)
    
%     %rename the arg list
%     num_args=length(func(i).arg_names);
%     arglist=[];
%     for k=1:num_args
%         func(i).tokenType(strcmpi(func(i).formulaToken,func(i).arg_names(k)))=7;
%         func(i).arg_names{k}=['arg',num2str(k)];
%         arglist=[arglist,'arg',num2str(k),','];
%     end
%     arglist=arglist(1:end-1); %remove last comma
%     
%     func(i).full_arglist=arglist;
    
    tokens=func(i).formulaToken(:);
    tokenType=func(i).tokenType;
    tokenIx=func(i).tokenIx;
%     thisFormula=buildFormula(tokens,tokenType,tokenIx,func(i));
%     output_file{end+1}=[name '=@(' func(i).full_arglist ') '  thisFormula ';'];
    output_file{end+1}=[name '=@(' func(i).full_arglist ') '   func(i).formula ';'];
end

for i=1:xppdata.nFixed
    
    name=fixed(i).name;
    
    tokens=fixed(i).formulaToken(:);
    tokenType=fixed(i).tokenType;
    tokenIx=fixed(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=[name '=' thisFormula ';'];
end

for i=1:xppdata.nVar
    
    if isempty(var(i).value), var(i).value=0;end
    
    tokens=var(i).formulaToken(:);
    tokenType=var(i).tokenType;
    tokenIx=var(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=[slopeName '(' num2str(i) ',:)=' thisFormula ';'];
end

for i=1:xppdata.nAux
    
    tokens=aux(i).formulaToken(:);
    tokenType=aux(i).tokenType;
    tokenIx=aux(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=[auxName '(' num2str(i) ',:)=' thisFormula ';'];
end

%Nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function thisFormula=buildFormula(tokens, tokenType, tokenIx, thisfunc)
        %Build formula in matlab syntax from tokens provided. Check the formula to
        %make any necessary changes
        
        
        thisFormula='';
        for j=1:length(tokens)
            
            tok=tokens{j};
            if tokenType(j)==0 %number
                
                %if a named number, get its value.
                if ischar(tok)
                    tok=num(tokenIx(j)).value;
                end
                
                %remove trailing zeros
                tok=sprintf('%g',tok);                
%                 tok=tok(1:find(tok~='0',1,'last'));
%                 tok=[tok,'0'];
                tok=['(' tok ')'];
                
            elseif tokenType(j)==1 %simple math
                %option to use vectorized ops 
                % - also needs convention for column/row-based ordering of
                % variables when provided as input (see tokentype=5)
                if options.vectorized
                    switch tok
                        case {'*','/','^'}
                            tok=['.',tok];
                    end
                end
                
            elseif tokenType(j)==2 %reserved math words. Convert to matlab version
                switch tok
                    
                    case 'flr'
                        tok='floor';
                        
                    case 'ln'
                        tok='log';
                        
                    case 'heav'
                        tok='heaviside';
                        
                    case 'lgamma'
                        tok='gammaln';
                        
                    case 'if'
                        error('if(condition)then(expression)else(expression)  is not yet implemented')
                        
                    case 'ran' %this takes more work. XPP: ran(maxNum) --> maxNum*rand()
                        error('ran is not yet implemented')
                        
                    case 'normal' %this takes more work. XPP: normal(mean,var)--> mean+sqrt(var)*randn()
                        error('normal is not yet implemented')
                        
                    case 'shift'
                        error('shift is not yet implemented')
                        
                    case 'del_shift'
                        error('del_shift is not yet implemented')
                        
                    case 'sum'
                        error('sum is not yet implemented')
                        
                    case 'hom_bcs'
                        error('hom_bcs is not yet implemented')
                end
                
                
            elseif tokenType(j)==3 %parameter
                tok=[paramName '(' num2str(tokenIx(j)) ')'];
                
                
            elseif tokenType(j)==4 %fixed quantity
                
                %nothing to do
                
            elseif tokenType(j)==5 %variable
%                 if options.vectorized
                    %convention: require row-wise vars like matlab odes
                    tok=[stateName '(' num2str(tokenIx(j)) ',:)'];
%                 else
%                     tok=[stateName '(' num2str(tokenIx(j)) ')'];
%                 end
                
            elseif tokenType(j)==6 %function name
                
                %nothing to change
                
            elseif tokenType(j)==7 %formal argument for a function
                
                %replace with new arg name
                tok=thisfunc.arg_names{tokenIx(j)};
                
            elseif tokenType(j)==8 %Wiener variable
                tok=[wienerName '(' num2str(tokenIx(j)) ')'];
                
            else
                error('Unidentified token type. Problem parsing??')
            end
            
            %accumulate the tokens into a single string
            thisFormula=[thisFormula, tok];
        end
    end

end