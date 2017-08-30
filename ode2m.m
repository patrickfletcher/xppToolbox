function [destinationFile, xppdata, clampVarIx]=ode2m(arg1, precision, destinationFile, clampVar)
% code conversion from ODE file for XPP to a .m file that can be run with
% MatLab. Calls parseODEfile, writes the ODE information to a .m file,
% then packages the secondary parameters into a vector pars. Also produced
% is a structure containing options found in the ODE file (e.g. dt, total,
% clODE stuff, etc.)
%
% USAGE
%   [mFunctionName, pars, y0, ivp, opt, clampVarIx]=ode2m(arg1, destinationFile, clampVar)
%
% Inputs
%   none - will prompt for ODE file, write to same name .m in working dir
%   arg1: either full path to ODE file, or pre-generated ivpStruct
%   destinationFile: specified output filename, NO extension (will be .m)
%   clampVar: name of variable to set as clamped variable
%
% Outputs
%   mFileName: string conatining the destination filename without extension
%   pars: default parameter values set in the ODE file
%   y0: default initial conditions set in the ODE file
%   opt: options structure contining all options found in ODE file (will be
%        empty if the ivpstruct was passed in as arg1)
%   clampVarIx: index of the variable that is clamped to its y0 value
%
% Example for using MatLab's ODE solvers:
%   [mFunctionName,p,y0]=ode2m('./xppSrc/lactotroph.ode');
%   fun=eval(['@(t,y)' mFunctionName '(t,y,p)']);
%   ode45(fun,[0,1000],y0);

%TO DO
% - standardize file I/O stuff!!!
% - both ode2m and ode2cl share several aspects... helper functions for input checking, I/O setup, etc??
% - destination file overwriting behavior? now destroys any existing copy
% - don't use symbolic toolbox heaviside!!

fileExtension='m';

% input checking
if nargin==0 || isempty(arg1)
    [name,path]=uigetfile('.ode','Select an ODE file');
    if ~ischar(name)
        disp('File selection canceled, quitting...')
        return
    end
    arg1=fullfile(path,name);
end

if ischar(arg1) %filename input: parse it.
    srcfilename=arg1;
    % Extract info from the ODE file
    xppdata=parseODEfile(srcfilename);
else %assume it is ivpStruct (TODO: errorchecking)
    xppdata=arg1;
end

%Build the destination filename
if ~exist('destinationFile','var')||isempty(destinationFile)
    destinationFile=xppdata.name;
else
    %destinationFile was provided. Extract just the filename without
    %extension or path.
    [~,destinationFile]=fileparts(destinationFile);
    
    %check whether it exists, and ask if overwriting it is ok?
end


if ~exist('precision','var') || isempty(precision)
    precision='double';
end


doClamp=false;
clampVarIx=[];
if exist('clampVar','var')&&~isempty(clampVar)
    %Check whether clamped variable name specified was in the ODE file
    clampVarIx=find(strcmpi({xppdata.var(:).name},clampVar));
    if ~isempty(clampVarIx)
        doClamp=true;
        destinationFile=[destinationFile '_' xppdata.var(clampVarIx).name '_clamp'];
    else
        disp('Warning: the variable specified for clamping was not found! No variables were clamped.')
    end
end

output_file=BuildOutputFile(xppdata,doClamp,precision,destinationFile);
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

%give suggestions on which matlab solver to try based on the method option
%detected. TODO: expand this with a better understanding of exactly which
%solvers XPP is using!

% if ~isempty(XPPopt.method)
%     switch XPPopt.method
%         case {'runge','rungekutta','modeuler','5dp','83dp'}
%             disp('Try ode45')
%         case {'backeuler','gear','stiff','cvode'}
%             disp('Try ode15s')
%         case '2rb'
%             disp('Try ode23s')
%     end
% end

end

function output_file=BuildOutputFile(xppdata,doClamp,precision,destinationFile)

if strcmpi(precision,'single')
    numFormat='%0.7f'; %7 digits, single precision floating point
elseif strcmpi(precision,'double')
    numFormat='%0.14f'; %14 digits, double precision floating point
end

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
stateName='state_vector';
paramName='param_vector';
slopeName='slope_vector';
auxName='aux_vector';

output_file={};
output_file{1}=['function [' slopeName ', ' auxName ']=' destinationFile '(t,' stateName ',' paramName ')'];
output_file{2}=[slopeName '=zeros(' num2str(xppdata.nVar) ',1,''' precision ''');'];
output_file{3}=[auxName '=zeros(' num2str(xppdata.nAux) ',1,''' precision ''');']; %initialize to empty if no Aux

%put anonymous function definitions first
for i=1:xppdata.nFunc
    
    name=func(i).name; %case is all lower after parser (XPP is case insensitive)
    
    %rename the arg list
    num_args=length(func(i).arg_names);
    arglist=[];
    for k=1:num_args
        func(i).tokenType(strcmpi(func(i).formulaToken,func(i).arg_names(k)))=7;
        func(i).arg_names{k}=['arg',num2str(k)];
        arglist=[arglist,'arg',num2str(k),','];
    end
    arglist=arglist(1:end-1); %remove last comma
    
    func(i).full_arglist=arglist;
    
    tokens=func(i).formulaToken(:);
    tokenType=func(i).tokenType;
    tokenIx=func(i).tokenIx;
    thisFormula=buildFormula(tokens,tokenType,tokenIx,func(i));
    output_file{end+1}=[name '=@(' func(i).full_arglist ') '  thisFormula ';'];
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
    
    if doClamp && i==clampVarIx
        output_file{end+1}=[slopeName '(' num2str(i) ',1)=0;'];
    else
        output_file{end+1}=[slopeName '(' num2str(i) ',1)=' thisFormula ';'];
    end
end

for i=1:xppdata.nAux
    
    tokens=aux(i).formulaToken(:);
    tokenType=aux(i).tokenType;
    tokenIx=aux(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=[auxName '(' num2str(i) ',1)=' thisFormula ';'];
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
                
                %specify format
                tok=sprintf(numFormat,tok);
                %remove trailing zeros
                tok=tok(1:find(tok~='0',1,'last'));
                tok=[tok,'0'];
                tok=['(' tok ')'];
                
            elseif tokenType(j)==1 %simple math
                
                %do nothing.
                
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
                tok=[stateName '(' num2str(tokenIx(j)) ')'];
                
            elseif tokenType(j)==6 %function name
                
                %nothing to change
                
            elseif tokenType(j)==7 %formal argument for a function
                
                %replace with new arg name
                tok=thisfunc.arg_names{tokenIx(j)};
            else
                error('Unidentified token type. Problem parsing??')
            end
            
            %accumulate the tokens into a single string
            thisFormula=[thisFormula, tok];
        end
    end

end