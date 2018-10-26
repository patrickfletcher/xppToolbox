function [destinationFile, xppdata]=xppConverter(sourceFile, destinationType, destinationFile, cSinglePrecision)
% code conversion from ODE file for XPP to other formats. These include:
% - Matlab: M or MEX files, for use in Matlab ODE solvers
% - clODE: CL file for use in clODE parallel ode solver
% - VFGEN: an XML format that can be used in VFGEN
% - VCML: an XML format compatible with VCML's math models.
%
% USAGE
%   [destinationFile, xppdata]=xppConverter(sourceFile, destinationType, destinationFile, precision)
%
% Inputs
%   none - will prompt for ODE file, write to same name .m in working dir
%   sourceFile - ode filename
%   destinationFile - output filename
%   precision - single/double, for c code
%
% Outputs
%   mFileName: string conatining the destination filename without extension
%   xppdata: structure containing all ODE file information
%
% Example for using MatLab's ODE solvers:
%   [mFunctionName,p,y0]=ode2m('./xppSrc/lactotroph.ode');
%   fun=eval(['@(t,y)' mFunctionName '(t,y,p)']);
%   ode45(fun,[0,1000],y0);

%TODO
% - refactor buildOutputFile subroutines: extract common functionality as
% single subroutines that can be called by any of them instead of nested
% subroutines for each.
% - update help blurb
% - destination file overwriting behavior? now destroys any existing copy
% - MEX: single precision doesn't work
% - MEX: package the template directly into this file?

% input checking
if nargin==0 || isempty(sourceFile)
    [name,path]=uigetfile('.ode','Select an ODE file');
    if ~ischar(name)
        disp('File selection canceled, quitting...')
        return
    end
    sourceFile=fullfile(path,name);
end

if ischar(sourceFile) %filename input: parse it.
    % Extract info from the ODE file
    xppdata=parseODEfile(sourceFile);
else %assume it is xppdata struct from previous parsing (TODO: errorchecking)
    xppdata=sourceFile;
end

%default: create a Matlab m-file
if ~exist('destinationType','var')||isempty(destinationType)
    destinationType='m';
end

%Get the destination file name (no extension yet)
if ~exist('destinationFile','var')||isempty(destinationFile)
    destinationFile=xppdata.name;
else
    %destinationFile was provided. Extract just the filename without
    %extension or path.
    [~,destinationFile]=fileparts(destinationFile);
    
    %check whether it exists, and ask if overwriting it is ok?
end

if ~exist('precision','var')
    cSinglePrecision=false;
end

switch lower(destinationType)
    case {'m'}
        fileExtension='m';
        output_file=BuildOutputFileM(xppdata,destinationFile);
        
    case {'mex'}
        fileExtension='c';
        output_file=BuildOutputFileMEX(xppdata,cSinglePrecision);
        
    case {'cl'}
        fileExtension='cl';
        output_file=BuildOutputFileCL(xppdata,cSinglePrecision);
        
    case {'vfgen'}
        fileExtension='vf';
        output_file=BuildOutputFileVFGEN(xppdata);
        
    case {'vcml'}
        fileExtension='vcml';
        output_file=BuildOutputFileVCML(xppdata);
        
    otherwise
        error('unrecognized destination type')
end

lineCount=length(output_file);

%delete a previous version - silently destroys any previous version!
fullDestinationPath=[pwd filesep destinationFile '.' fileExtension];
if exist(fullDestinationPath,'file')==2
    delete(fullDestinationPath)
end

% write to file
fidw=fopen(fullDestinationPath,'w');

if fidw > -1
    for line=1:lineCount
        fprintf(fidw, '%s\n', output_file{line});
    end
    status=fclose(fidw);
else
    fprintf(' Problem opening file %s for writing. Cannot continue\n',fullDestinationPath);
    return
end

%check to make sure file is written
ready=false;
while ~ready
    status=exist(fullDestinationPath,'file');
    if status
        ready=true;
    end
end


% type-specific cleanup
switch lower(destinationType)
    case {'m'}
        
    case {'mex'}
        %compile the mex file
        mex(fullDestinationPath);
        disp('mex file generated')
        
    case {'cl'}
        xppdata.clRHSfilename=fullDestinationPath;
        destinationFile=fullDestinationPath;
        
    case {'vfgen'}
        
    case {'vcml'}
        
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


%%% Helper functions %%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BuildOutputFileM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_file=BuildOutputFileM(xppdata,destinationFile)

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
output_file{3}=[auxName '=zeros(' num2str(xppdata.nAux) ',1);']; %initialize to empty if no Aux

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
                
                %specify format
                tok=sprintf('%g',tok);
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


%%% C-language helpers
function [tokens,tokenType,tokenIx]=convertPowers(tokens,tokenType,tokenIx)
%convert "a^b" to "pow(a,b)" - find a(:),^,b(:) then replace with
%{'pow' '(' a(:) ',' b(:) ')'}

%TODO: find a way to check if exponent is integer, then use pown
%if isreal(x) && rem(x,1)==0, where x=expon{:} concatenated... If only
%numbers are in expon, should work, but if names, would need to recursively
%check if it all boils down to a number

nPow=sum(strcmpi(tokens,'^'));
for j=1:nPow
    base={};baseType=[];baseTokIx=[];
    expon={};exponType=[];exponTokIx=[];
    thisPowIx=find(strcmpi(tokens,'^'),1,'first');
    %find base: number/name or something in parentheses
    if tokenType(thisPowIx-1)==1
        %assume it is ')' for now: full validity of math expressions is not checked!!!!
        base={tokens{thisPowIx-1}};
        baseType=tokenType(thisPowIx-1);
        baseTokIx=tokenIx(thisPowIx-1);
        
        parenCount=1; %increment for ')', decrement for '('
        ix=1;
        while parenCount>0
            ix=ix+1;
            %check whether parentheses have closed
            if strcmp(tokens(thisPowIx-ix),'(')
                parenCount=parenCount-1;
            elseif strcmp(tokens(thisPowIx-ix),')')
                parenCount=parenCount+1;
            end
            
            %add the new token to beginning of base
            base={tokens{thisPowIx-ix},base{:}};
            baseType=[tokenType(thisPowIx-ix),baseType];
            baseTokIx=[tokenIx(thisPowIx-ix),baseTokIx];
            baseStartIx=thisPowIx-ix;
        end
        
    else
        base={tokens{thisPowIx-1}};
        baseType=tokenType(thisPowIx-1);
        baseTokIx=tokenIx(thisPowIx-1);
        baseStartIx=thisPowIx-1;
    end
    
    %find exponent: number/name or something in parentheses
    if tokenType(thisPowIx+1)==1
        
        expon={tokens{thisPowIx+1}};
        exponType=tokenType(thisPowIx+1);
        exponTokIx=tokenIx(thisPowIx+1);
        
        parenCount=1; %increment for '(', decrement for ')'
        ix=1;
        while parenCount>0
            ix=ix+1;
            %check whether parentheses have closed
            if strcmp(tokens(thisPowIx+ix),')')
                parenCount=parenCount-1;
            elseif strcmp(tokens(thisPowIx+ix),'(')
                parenCount=parenCount+1;
            end
            
            %add the new token to beginning of base
            expon={expon{:},tokens{thisPowIx+ix}};
            exponType=[exponType,tokenType(thisPowIx+ix)];
            exponTokIx=[exponTokIx,tokenIx(thisPowIx+ix)];
            exponEndIx=thisPowIx+ix;
        end
        
    else
        expon={tokens{thisPowIx+1}};
        exponType=tokenType(thisPowIx+1);
        exponTokIx=tokenIx(thisPowIx+1);
        exponEndIx=thisPowIx+1;
    end
    
    tokens={tokens{1:baseStartIx-1},'pow','(', base{:}, ',', expon{:}, ')', tokens{exponEndIx+1:end} };
    tokenType=[tokenType(1:baseStartIx-1),1,1,baseType,1,exponType,1,tokenType(exponEndIx+1:end)];
    tokenIx=[tokenIx(1:baseStartIx-1),0,0,baseTokIx,0,exponTokIx,0,tokenIx(exponEndIx+1:end)];
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BuildOutputFileMEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_file=BuildOutputFileMEX(xppdata,cSinglePrecision)

if cSinglePrecision
    numFormat='%0.7g'; %7 digits, single precision floating point TODO: check correctness
else
    numFormat='%0.14g'; %14 digits, double precision floating point
end

if xppdata.nWiener>0
    error('Wiener process not yet supported in ode2mex')
end

num=xppdata.num;
var=xppdata.var;
fixed=xppdata.fixed;
func=xppdata.func;
aux=xppdata.aux;

%names for internal vectors, to try to avoid potential name collisions
%with user fixed quantities, parameters, or variables.
stateName='x_';
paramName='p_';
slopeName='dx_';
auxName='aux_';

% build the output as a cell array containing each line.
% Structure:
%   1. Defines (for error checking mex inputs, controlling precision)
%   2. Inline function definitions
%   3. function signature
%   4. initialize output variables
%   5. fixed formulas
%   6. differential equations
%   7. auxiliary output variables

output_file{1}='#include <math.h>';
output_file{end+1}='#include "mex.h"';
if cSinglePrecision
    output_file{end+1}='#define realtype float';
else
    output_file{end+1}='#define realtype double'; %14 digits, double precision floating point
end

output_file{end+1}=['#define N_PAR ' num2str(xppdata.nPar)];
output_file{end+1}=['#define N_VAR ' num2str(xppdata.nVar)];
output_file{end+1}=['#define N_AUX ' num2str(max(xppdata.nAux,1))]; %number of auxiliary variables should be at least 1... bug?

%function notation not yet supported for OpenCL.
if xppdata.nFunc>0
    disp('Error processing line: ')
    disp([func(1).name '(' func(1).full_arglist ')='  func(1).formula])
    error('"name(var1,var2,...)=formula" function notation not yet supported');
end
% %put function definitions first
% for i=1:xppdata.nFunc
%     
%     name=func(i).name; %case is all lower after parser (XPP is case insensitive)
%     
%     %rename the arg list
%     %     num_args=length(func(i).arg_names);
%     %     arglist=[];
%     %     for k=1:num_args
%     %         func(i).arg_names{k}=['arg',num2str(k)];
%     %         arglist=[arglist,'realtype arg',num2str(k),','];
%     %     end
%     %     arglist=arglist(1:end-1); %remove last comma
%     %
%     %     func(i).full_arglist=arglist;
%     
%     tokens=func(i).formulaToken(:);
%     tokenType=func(i).tokenType;
%     tokenIx=func(i).tokenIx;
%     thisFormula=buildFormula(tokens,tokenType,tokenIx,func(i));
%     
%     %     output_file{end+1}=['#define ' name '(' func(i).full_arglist ') ' thisFormula];
%     
%     %     output_file{end+1}=['realtype ' name '(realtype t, realtype p[], realtype y[], ' func(i).full_arglist ') {' ];
%     %     output_file{end+1}=['   return ' thisFormula ';'];
%     %     output_file{end+1}=['}'];
%     output_file{end+1}='';
% end

%main body of RHS file
output_file{end+1}=['static void yprime(realtype *t, realtype ' ...
    stateName '[], realtype ' paramName '[], realtype ' slopeName '[], realtype ' auxName '[]) {'];

for i=1:xppdata.nFixed
    
    name=fixed(i).name;
    
    tokens=fixed(i).formulaToken(:);
    tokenType=fixed(i).tokenType;
    tokenIx=fixed(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=['realtype ' name '=' thisFormula ';'];
    
end

for i=1:xppdata.nVar
    
    tokens=var(i).formulaToken(:);
    tokenType=var(i).tokenType;
    tokenIx=var(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    output_file{end+1}=[slopeName '[' num2str(i-1) ']=' thisFormula ';'];
end

for i=1:xppdata.nAux
    
    tokens=aux(i).formulaToken(:);
    tokenType=aux(i).tokenType;
    tokenIx=aux(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=[auxName '[' num2str(i-1) ']=' thisFormula ';'];
end

output_file{end+1}='}';
output_file{end+1}='';

% now collect the lines from the mexRHStemplate
fidt=fopen('mexRHStemplate.c','r');
tline = fgetl(fidt);
while ischar(tline)
    output_file{end+1}=tline;
    tline=fgetl(fidt);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function thisFormula=buildFormula(tokens, tokenType, tokenIx, thisfunc)
        
        [tokens,tokenType,tokenIx]=convertPowers(tokens,tokenType,tokenIx);
        
        %functions of form: name(arg1,arg2,...,arg3) can be left as is, and the
        %args modified as appropriate, as long as "name" is a known predefined
        %helper function like heaviside, mod
        
        %check the formula to make any necessary changes
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
%                 tok=tok(1:find(tok~='0',1,'last'));
%                 tok=[tok,'0'];
                
                if cSinglePrecision %add an f to indicate single precision
                    tok=[tok,'f'];
                end
                
                tok=['(' tok ')'];
                
            elseif tokenType(j)==1 %simple math
                
                %nothing to do.
                
            elseif tokenType(j)==2 %reserved math word
                
                switch tok
                    case 'abs'
                        tok='fabs';
                        
                        %             case 'heav' %helper function named heav is defined
                        %             tok='step';
                        
                        
                    case 't'
                        tok='*t';
                        
                    case 'flr'
                        tok='floor';
                        
                    case 'ln'
                        tok='log';
                        
                    case 'max'
                        tok='fmax';
                        
                    case 'min'
                        tok='fmin';
                        
                    case 'mod'
                        tok='fmod';
                        
                    case 'pi'
                        if strcmpi(precision,'single')
                            tok='M_PI_F';
                        elseif strcmpi(precision,'double')
                            tok='M_PI';
                        end
                        
                        
                    case {'besseli', 'besselj', 'bessely'}
                        %                 tok='rand()';
                        error('Bessel functions not yet implemented')
                        
                    case 'if'
                        %                 tok='rand()';
                        error('if(condition)then(expression)else(expression)  is not yet implemented')
                        
                    case 'ran' %this takes more work. XPP: ran(maxNum) --> maxNum*rand()
                        %                 tok='rand()';
                        error('ran is not yet implemented')
                        
                    case 'normal' %this takes more work. XPP: normal(mean,var)--> mean+sqrt(var)*randn()
                        %                 tok='randn()';
                        error('normal is not yet implemented')
                        
                    case 'shift'
                        %                 tok='rand()';
                        error('shift is not yet implemented')
                        
                    case 'sign'
                        %                 tok='rand()';
                        error('sign is not yet implemented')
                        
                    case 'del_shift'
                        %                 tok='rand()';
                        error('del_shift is not yet implemented')
                        
                    case 'sum'
                        %                 tok='rand()';
                        error('sum is not yet implemented')
                        
                    case 'hom_bcs'
                        %                 tok='rand()';
                        error('hom_bcs is not yet implemented')
                        
                end
                
            elseif tokenType(j)==3 %parameter
                tok=[paramName '[' num2str(tokenIx(j)-1) ']'];
                
            elseif tokenType(j)==4 %fixed quantity
                %do nothing
                
            elseif tokenType(j)==5 %variable
                %replace with 'y(k)' where k=variable index
                tok=[stateName '[' num2str(tokenIx(j)-1) ']'];
                
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BuildOutputFileCL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_file=BuildOutputFileCL(xppdata,cSinglePrecision)

if cSinglePrecision
    numFormat='%0.7g'; %7 digits, single precision floating point TODO: check correctness
else
    numFormat='%0.14g'; %14 digits, double precision floating point
end

num=xppdata.num;
var=xppdata.var;
fixed=xppdata.fixed;
func=xppdata.func;
aux=xppdata.aux;


% build the output as a cell array containing each line.
% Structure:
%   0. Inline function definitions
%   1. function signature
%   2. initialize output variables
%   3. fixed formulas
%   4. differential equations
%   5. auxiliary output variables


%names for internal vectors, to try to avoid potential name collisions
%with user fixed quantities, parameters, or variables.
stateName='x_';
paramName='p_';
slopeName='dx_';
auxName='aux_';
wienerName='w_';

%the opencl header defines realtype to be float or double, based on precision

output_file{1}='';

%function notation not yet supported for OpenCL.
if xppdata.nFunc>0
    disp('Error processing line: ')
    disp([func(1).name '(' func(1).full_arglist ')='  func(1).formula])
    error('"name(var1,var2,...)=formula" function notation not yet supported');
end

% % %put function definitions first
% % for i=1:xppdata.nFunc
% %
% %     name=func(i).name; %case is all lower after parser (XPP is case insensitive)
% %
% %     %rename the arg list
% % %     num_args=length(func(i).arg_names);
% % %     arglist=[];
% % %     for k=1:num_args
% % %         func(i).arg_names{k}=['arg',num2str(k)];
% % %         arglist=[arglist,'realtype arg',num2str(k),','];
% % %     end
% % %     arglist=arglist(1:end-1); %remove last comma
% % %
% % %     func(i).full_arglist=arglist;
% %
% %     tokens=func(i).formulaToken(:);
% %     tokenType=func(i).tokenType;
% %     tokenIx=func(i).tokenIx;
% %     thisFormula=buildFormula(tokens,tokenType,tokenIx,func(i));
% %
% % %     output_file{end+1}=['#define ' name '(' func(i).full_arglist ') ' thisFormula];
% %
% % %     output_file{end+1}=['realtype ' name '(realtype t, realtype p[], realtype y[], ' func(i).full_arglist ') {' ];
% % %     output_file{end+1}=['   return ' thisFormula ';'];
% % %     output_file{end+1}=['}'];
% %     output_file{1}='';
% % end

%main body of RHS file
output_file{end+1}=['void getRHS(realtype t, realtype ' ...
    stateName '[], realtype ' paramName '[], realtype ' slopeName '[], realtype ' auxName '[], realtype ' wienerName '[]) {'];

for i=1:xppdata.nFixed
    
    name=fixed(i).name;
    
    tokens=fixed(i).formulaToken(:);
    tokenType=fixed(i).tokenType;
    tokenIx=fixed(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=['realtype ' name '=' thisFormula ';'];
    
end

for i=1:xppdata.nVar
    
    tokens=var(i).formulaToken(:);
    tokenType=var(i).tokenType;
    tokenIx=var(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    output_file{end+1}=[slopeName '[' num2str(i-1) ']=' thisFormula ';'];
end

for i=1:xppdata.nAux
    
    tokens=aux(i).formulaToken(:);
    tokenType=aux(i).tokenType;
    tokenIx=aux(i).tokenIx;
    
    thisFormula=buildFormula(tokens,tokenType,tokenIx);
    
    output_file{end+1}=[auxName '[' num2str(i-1) ']=' thisFormula ';'];
end

output_file{end+1}='}';
output_file{end+1}='';


%Nested function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function thisFormula=buildFormula(tokens, tokenType, tokenIx, thisfunc)
        
        %TODO: discover powers with integer exponents, convert to intPower(a,b)
        
        [tokens,tokenType,tokenIx]=convertPowers(tokens,tokenType,tokenIx);
        
        %functions of form: name(arg1,arg2,...,arg3) can be left as is, and the
        %args modified as appropriate, as long as "name" is a known predefined
        %helper function like heaviside, mod
        
        %check the formula to make any necessary changes
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
%                 tok=tok(1:find(tok~='0',1,'last'));
%                 tok=[tok,'0'];
                
                if cSinglePrecision %add an f to indicate single precision
                    tok=[tok,'f'];
                end
                
                tok=['(' tok ')'];
                
            elseif tokenType(j)==1 %simple math
                
                %nothing to do.
                
            elseif tokenType(j)==2 %reserved math word
                
                switch tok
                    case 'abs'
                        tok='fabs';
                        
                        %             case 'heav' %helper function named heav is defined
                        %             tok='step';
                        
                        
                    case 'flr'
                        tok='floor';
                        
                    case 'ln'
                        tok='log';
                        
                    case 'max'
                        tok='fmax';
                        
                    case 'min'
                        tok='fmin';
                        
                    case 'mod'
                        tok='fmod';
                        
                    case 'pi'
                        if strcmpi(precision,'single')
                            tok='M_PI_F';
                        elseif strcmpi(precision,'double')
                            tok='M_PI';
                        end
                        
                        
                    case {'besseli', 'besselj', 'bessely'}
                        %                 tok='rand()';
                        error('Bessel functions not yet implemented')
                        
                    case 'if'
                        %                 tok='rand()';
                        error('if(condition)then(expression)else(expression)  is not yet implemented')
                        
                    case 'ran' %this takes more work. XPP: ran(maxNum) --> maxNum*rand()
                        %                 tok='rand()';
                        error('ran is not yet implemented')
                        
                    case 'normal' %this takes more work. XPP: normal(mean,var)--> mean+sqrt(var)*randn()
                        %                 tok='randn()';
                        error('normal is not yet implemented')
                        
                    case 'shift'
                        %                 tok='rand()';
                        error('shift is not yet implemented')
                        
                    case 'sign'
                        %                 tok='rand()';
                        error('sign is not yet implemented')
                        
                    case 'del_shift'
                        %                 tok='rand()';
                        error('del_shift is not yet implemented')
                        
                    case 'sum'
                        %                 tok='rand()';
                        error('sum is not yet implemented')
                        
                    case 'hom_bcs'
                        %                 tok='rand()';
                        error('hom_bcs is not yet implemented')
                        
                end
                
            elseif tokenType(j)==3 %parameter
                tok=[paramName '[' num2str(tokenIx(j)-1) ']'];
                
            elseif tokenType(j)==4 %fixed quantity
                %do nothing
                
            elseif tokenType(j)==5 %variable
                tok=[stateName '[' num2str(tokenIx(j)-1) ']'];
                
            elseif tokenType(j)==6 %function name
                %nothing to change
                
            elseif tokenType(j)==7 %formal argument for a function
                %replace with new arg name
                tok=thisfunc.arg_names{tokenIx(j)};
                
            elseif tokenType(j)==8 %Wiener variable
                tok=[wienerName '[' num2str(tokenIx(j)-1) ']'];
                
            else
                error('Unidentified token type. Problem parsing??')
            end
            
            %accumulate the tokens into a single string
            thisFormula=[thisFormula, tok];
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BuildOutputFileVFGEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_file=BuildOutputFileVFGEN(xppdata)

% build the output as a cell array containing each line.
par=xppdata.par;
num=xppdata.num;
var=xppdata.var;
fixed=xppdata.fixed;
func=xppdata.func;
aux=xppdata.aux;

if xppdata.nFunc>0
    disp([func(1).name '(' func(1).full_arglist ')='  func(1).formula])
    error('"name(var1,var2,...)=formula" function notation not yet supported');
end


output_file={};
output_file{end+1}='<?xml version="1.0" ?>';
output_file{end+1}='<!--This file was generated by ode2vfgen-->';


%Begin VectorField block
output_file{end+1}=['<VectorField Name="' xppdata.name '">'];

%constants
for i=1:xppdata.nNum
    output_file{end+1}=['<Constant Name="' num(i).name '" Value="' num2str(num(i).value) '" />'];
end

%parameters
for i=1:xppdata.nPar
    output_file{end+1}=['<Parameter Name="' par(i).name '" DefaultValue="' num2str(par(i).value) '" />'];
end

%Expressions
for i=1:xppdata.nFixed
    output_file{end+1}=['<Expression  Name="' fixed(i).name '" Formula="' fixed(i).formula '"/>'];
end

%State Variables
for i=1:xppdata.nVar
    output_file{end+1}=['<StateVariable   Name="' var(i).name '" ' ...
                                         'Formula="' var(i).formula '" '...
                                         'DefaultInitialCondition="' num2str(var(i).value) '"/>'];
end

%User Functions
for i=1:xppdata.nAux
    output_file{end+1}=['<Function  Name="aux_' aux(i).name '" Formula="' aux(i).formula '"/>'];
end

output_file{end+1}='</VectorField>';
%End VectorField Block
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BuildOutputFileVCML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_file=BuildOutputFileVCML(xppdata)

par=xppdata.par;
num=xppdata.num;
var=xppdata.var;
fixed=xppdata.fixed;
func=xppdata.func;
aux=xppdata.aux;
opt=xppdata.opt;

if xppdata.nFunc>0
    disp([func(1).name '(' func(1).full_arglist ')='  func(1).formula])
    error('"name(var1,var2,...)=formula" function notation not yet supported');
end


% build the output as a cell array, one line per element
output_file={};
output_file{end+1}='<?xml version="1.0" encoding="UTF-8"?>';
output_file{end+1}='<!--This file was generated by ode2vcml-->';

%Begin VCML block
output_file{end+1}='<vcml xmlns="http://sourceforge.net/projects/vcell/vcml" Version="Rel_Version_5.3_build_14">'; %

%Begin MathModel block
output_file{end+1}=['<MathModel Name="' xppdata.name '">'];

output_file{end+1}=['<Annotation>' xppdata.comment(1).text '</Annotation>'];

%dummy geometry block
output_file{end+1}='<Geometry Name="noname" Dimension="0">';
output_file{end+1}='<Extent X="10.0" Y="10.0" Z="10.0"/>';
output_file{end+1}='<Origin X="0.0" Y="0.0" Z="0.0"/>';
output_file{end+1}='<SubVolume Name="Compartment" Handle="0" Type="Compartmental"/>';
output_file{end+1}='</Geometry>';

%Begin MathDescription Block
output_file{end+1}='<MathDescription Name="noname">';

%constants
for i=1:xppdata.nPar
    checkVCMLreservedNames(par(i).name)
    output_file{end+1}=['<Constant Name="' par(i).name '">' num2str(par(i).value) '</Constant>'];
end

for i=1:xppdata.nNum
    checkVCMLreservedNames(num(i).name)
    output_file{end+1}=['<Constant Name="' num(i).name '">' num2str(num(i).value) '</Constant>'];
end

%generate new names for initial considions
for i=1:xppdata.nVar
    checkVCMLreservedNames(var(i).name)
    output_file{end+1}=['<Constant  Name="init_' var(i).name '">' num2str(var(i).value) '</Constant>'];
end

%volume varible names
for i=1:xppdata.nVar
    output_file{end+1}=['<VolumeVariable  Name="' var(i).name '"/>'];
end

%Functions - convert to fixed? need to find all instances of USE of
%function, and make a fixed expression for each...

% for i=1:nFunc
%     output_file{end+1}=['<Function  Name="' func(i).name '">' func(i).formula '</Function>'];
% end

for i=1:xppdata.nFixed
    checkVCMLreservedNames(fixed(i).name)
    fixed(i).formula=checkFormula(fixed(i).formula);
    output_file{end+1}=['<Function  Name="' fixed(i).name '">' fixed(i).formula '</Function>'];
end

for i=1:xppdata.nAux
    checkVCMLreservedNames(aux(i).name)
    aux(i).formula=checkFormula(aux(i).formula);
    output_file{end+1}=['<Function  Name="aux_' aux(i).name '">' aux(i).formula '</Function>'];
end


%Begin Compartment Subdomain
output_file{end+1}='<CompartmentSubDomain Name="Compartment">';
output_file{end+1}='<BoundaryType Boundary="Xm" Type="Value"/>';
output_file{end+1}='<BoundaryType Boundary="Xp" Type="Value"/>';
output_file{end+1}='<BoundaryType Boundary="Ym" Type="Value"/>';
output_file{end+1}='<BoundaryType Boundary="Yp" Type="Value"/>';
output_file{end+1}='<BoundaryType Boundary="Zm" Type="Value"/>';
output_file{end+1}='<BoundaryType Boundary="Zp" Type="Value"/>';

%ODE equations
for i=1:xppdata.nVar
    var(i).formula=checkFormula(var(i).formula);
    output_file{end+1}=['<OdeEquation   Name="' var(i).name '" SolutionType="Unknown">'];
    output_file{end+1}=['<Rate>' var(i).formula '</Rate>'];
    output_file{end+1}=['<Initial>init_' var(i).name '</Initial></OdeEquation>'];
end


output_file{end+1}='</CompartmentSubDomain>';
%End Compartment Subdomain


output_file{end+1}='</MathDescription>';
%End MathDescription Block


%Begin Simulation Block

%extract solvername from ode file. Try to map it to one of the following
%options:
% Forward Euler,
% Runge-Kutta (Second Order, Fixed Time Step),
% Runge-Kutta (Fourth Order, Fixed Time Step),
% Runge-Kutta-Fehlberg (Fifth Order, Variable Time Step),
% Adams-Moulton (Fifth Order, Fixed Time Step),
% IDA (Variable Order, Variable Time Step, ODE/DAE),
% CVODE (Variable Order, Variable Time Step).
%
% Default: Runge-Kutta-Fehlberg (Fifth Order, Variable Time Step)

%need to work out what aliases XPP allows for each, and if it is even
%supported

switch opt.method
    case 'euler'
        solvername='Forward Euler';
    case 'modeuler'
        solvername='Runge-Kutta (Second Order, Fixed Time Step)';
    case 'runge'
        solvername='Runge-Kutta (Fourth Order, Fixed Time Step)';
    case 'qualrk'
        solvername='Runge-Kutta-Fehlberg (Fifth Order, Variable Time Step)';
    case 'adams'
        solvername='Adams-Moulton (Fifth Order, Fixed Time Step)';
        %     case 'ida'
        %         name='IDA (Variable Order, Variable Time Step, ODE/DAE)';
    case 'cvode'
        solvername='CVODE (Variable Order, Variable Time Step)';
    otherwise
        solvername='Runge-Kutta-Fehlberg (Fifth Order, Variable Time Step)';
end

output_file{end+1}='<Simulation Name="Simulation1">';

%note: these will be populated by defaults set in parseODEfile if not set
%in the ode file directly
output_file{end+1}=['<SolverTaskDescription TaskType="Unsteady" UseSymbolicJacobian="false" Solver="' solvername '">'];
output_file{end+1}=['<TimeBound StartTime="' num2str(opt.t0) '" EndTime="' num2str(opt.total) '"/>'];
output_file{end+1}=['<TimeStep DefaultTime="' num2str(opt.dt) '" MinTime="' num2str(opt.dtmin) '" MaxTime="' num2str(opt.dtmax) '" />'];
output_file{end+1}=['<ErrorTolerance Absolut="' num2str(opt.atoler) '" Relative="' num2str(opt.toler) '" />'];
output_file{end+1}=['<OutputOptions KeepEvery="' num2str(opt.nout) '" KeepAtMost="' num2str(opt.maxstor) '" />'];
output_file{end+1}='<NumberProcessors>1</NumberProcessors>';
output_file{end+1}='</SolverTaskDescription>';

output_file{end+1}='<MathOverrides/>';

%For whatever reason, the simulation needs to have version info
output_file{end+1}='<Version Name="Simulation1" KeyValue="99999999" BranchId="99999999" Archived="0" Date="1-Jan-2312 00:00:00" FromVersionable="false">';
output_file{end+1}='<Owner Name="" Identifier="99999999"/>';
output_file{end+1}='<GroupAccess Type="1"/>';
output_file{end+1}='<Annotation> generated by ode2vcml </Annotation>';
output_file{end+1}='</Version>';


output_file{end+1}='</Simulation>';
% %End Simulation Block


output_file{end+1}='</MathModel>';
%End MathDescription Block


output_file{end+1}='</vcml>';
%End MathDescription Block


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Helper functions to convert to VCML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function checkVCMLreservedNames(name)
        reserved={'x','y','z','t'};
        if ismember(name,reserved)
            disp('')
            disp(['---> ' name])
            error('The above name is reserved by Virtual Cell. Please rename it and try again.')
        end
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function formula=checkFormula(formula)
        %Find unsupported RHS elements.
        % throw errors for unsupported special xpp names
        
        specialXPPnames = {'heav(', 'sign(', 'mod(',...
            'delay(', 'del_shft(', 'shift(', 'ishift(',...
            'ran(', 'normal(','lgamma(', ...
            'besseli(', 'besselj(', 'bessely(', 'erf(', 'erfc(', ...
            'hom_bcs(',...
            };
        
        %This is a bit of a bandaid right now, would be better to tokenize and
        %classify each element of the formula. This would allow more robust
        %replacement of specific tokens (e.g. pi-->3.1415...)
        % 'ln' ->log?
        % 'flr' -> floor
        % 'pi' -> numerical value
        
        if strfind(formula,'ln(')
            formula=strrep(formula,'ln(','log(');
        end
        if strfind(formula,'flr(')
            formula=strrep(formula,'flr(','floor(');
        end
        
        
        for j=1:length(specialXPPnames)
            if strfind(formula,specialXPPnames{j})
                disp(['--> ' formula])
                error([specialXPPnames{j}(1:end-1) ' not supported'])
            end
        end
    end


end