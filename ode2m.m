function [destinationFile, xppdata]=ode2m(source, destinationFile)
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
%   from parser
%   destinationFile: specified output filename
%
% Outputs
%   mFileName: string conatining the destination filename without extension
%   xppdata: structure containing all ODE file information
%
% Example for using MatLab's ODE solvers:
%   [mFunctionName,p,y0]=ode2m('./xppSrc/lactotroph.ode');
%   fun=eval(['@(t,y)' mFunctionName '(t,y,p)']);
%   ode45(fun,[0,1000],y0);

%TO DO
% - destination file overwriting behavior? now destroys any existing copy
% - don't use symbolic toolbox heaviside!!

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

%Build the destination filename
if ~exist('destinationFile','var')||isempty(destinationFile)
    destinationFile=xppdata.name;
else
    %destinationFile was provided. Extract just the filename without
    %extension or path.
    [~,destinationFile]=fileparts(destinationFile);
    
    %check whether it exists, and ask if overwriting it is ok?
end

output_file=BuildOutputFile(xppdata,destinationFile);
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

function output_file=BuildOutputFile(xppdata,destinationFile)

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

output_file={};
if xppdata.nWiener==0
output_file{1}=['function [' slopeName ', ' auxName ']=' destinationFile '(t,' stateName ',' paramName ')'];
else
output_file{1}=['function [' slopeName ', ' auxName ']=' destinationFile '(t,' stateName ',' wienerName ',' paramName ')'];
end
output_file{2}=[slopeName '=zeros(' num2str(xppdata.nVar) ',1);'];
output_file{3}=[auxName '=zeros(' num2str(xppdata.nAux) ',1);'];

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
    
    output_file{end+1}=[slopeName '(' num2str(i) ',1)=' thisFormula ';'];
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
                
                %remove trailing zeros
                tok=sprintf('%g',tok);
                
%                 tok=tok(1:find(tok~='0',1,'last'));
%                 tok=[tok,'0'];
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