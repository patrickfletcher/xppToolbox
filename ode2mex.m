function [destinationFile, xppdata, clampVarIx]=ode2m(arg1, precision, destinationFile, clampVar)

% code conversion from ODE file for XPP to a MEX file that can be run with
% MatLab ODE solvers. Calls parseXPPodeFile, writes the ODE information to
% a .c file, then compiles this into a mex file in the current directory.
% The secondary parameters are packaged into a vector pars. Also produced
% is a structure containing options found in the ODE file (e.g. dt, total,
% clODE stuff, etc.)
%
% this will make a C function that computes the RHS for the system given in
% the ODE file.
%
% specifying the destination filename is optional; if not specified, the
% output filename will be same as input filename, but with .c suffix in
% the current folder
%
% pars output contains default par values
% y0 output contains initial variable values


%use ccoder(sym(expression)) for getting power/sqrt etc right?

fileExtension='c';
opt=struct('name',{},'value',{});

% input checking
if nargin==0
    [name,path]=uigetfile('.ode','Select an ODE file');
    arg1=fullfile(path,name);
end

if ischar(arg1) %filename input: parse it.
    srcfilename=arg1;
    % Extract info from the ODE file
    xppdata=parseODEfile(srcfilename);
else %assume it is ivpStruct (TODO: errorchecking)
    xppdata=arg1;
end

ivpname=xppdata.name;

%Build the destination filename
if ~exist('destinationFile','var')||isempty(destinationFile)
    destinationFile=ivpname;
else
    fileSepIx=find(destinationFile==filesep,1,'last');
    if isempty(fileSepIx)
        mFunctionName=destinationFile;
    else
        mFunctionName=destinationFile(fileSepIx+1:end);
    end
end

doClamp=false;
clampVarIx=[];
if exist('clampVar','var')&&~isempty(clampVar)
    %Check whether clamped variable name specified was in the ODE file
    clampVarIx=find(strcmpi({var(:).name},clampVar));
    if ~isempty(clampVarIx)
        doClamp=true;
        destinationFile=[destinationFile '_' var(clampVarIx).name '_clamp'];
    else
        disp('Warning: the variable specified for clamping was not found! No variables were clamped.')
    end
end

% destinationFile=[destinationFile '.' fileExtension];
destinationFile=[destinationFile '_mex.' fileExtension]; %for testing, so i can also have the m-file
%end build destination filename

if ~exist('precision','var') || isempty(precision)
    precision='double';
end


output_file=BuildOutputFile(xppdata,doClamp,precision,destinationFile);
lineCount=length(output_file);


%delete a previous version - silently destroys any previous version!
fullPath=[pwd filesep destinationFile];
if exist(fullPath,'file')==2
    delete(fullPath)
end

% write to file
fidw=fopen(fullPath,'w');

if fidw ~= -1
    for line=1:lineCount
        fprintf(fidw, '%s\r\n', output_file{line});
    end
    
    fclose(fidw);
else
    fprintf(' Problem opening file %s for writing. Cannot continue\n',fullPath);
    return
end

%compile the mex file
mex(fullPath);
disp('mex file generated')


[~,destinationFile]=fileparts(destinationFile);

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

%names for internal vectors, to try to avoid potential name collisions
%with user fixed quantities, parameters, or variables.
stateName='state_vector';
paramName='param_vector';
slopeName='slope_vector';
auxName='aux_vector';

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
if strcmpi(precision,'single')
    output_file{end+1}='#define realtype float';
elseif strcmpi(precision,'double')
    output_file{end+1}='#define realtype double'; %14 digits, double precision floating point
end

output_file{end+1}=['#define N_PAR ' num2str(xppdata.nPar)];
output_file{end+1}=['#define N_VAR ' num2str(xppdata.nVar)];
output_file{end+1}=['#define N_AUX ' num2str(max(xppdata.nAux,1))]; %number of auxiliary variables should be at least 1... bug?

%put function definitions first
for i=1:xppdata.nFunc
    
    name=func(i).name; %case is all lower after parser (XPP is case insensitive)
    
    %rename the arg list
    %     num_args=length(func(i).arg_names);
    %     arglist=[];
    %     for k=1:num_args
    %         func(i).arg_names{k}=['arg',num2str(k)];
    %         arglist=[arglist,'realtype arg',num2str(k),','];
    %     end
    %     arglist=arglist(1:end-1); %remove last comma
    %
    %     func(i).full_arglist=arglist;
    
    tokens=func(i).formulaToken(:);
    tokenType=func(i).tokenType;
    tokenIx=func(i).tokenIx;
    thisFormula=buildFormula(tokens,tokenType,tokenIx,func(i));
    
    %     output_file{end+1}=['#define ' name '(' func(i).full_arglist ') ' thisFormula];
    
    %     output_file{end+1}=['realtype ' name '(realtype t, realtype p[], realtype y[], ' func(i).full_arglist ') {' ];
    %     output_file{end+1}=['   return ' thisFormula ';'];
    %     output_file{end+1}=['}'];
    output_file{end+1}='';
end

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
    
    if doClamp && i==clampVarIx
        output_file{end+1}=[slopeName '[' num2str(i-1) ']=0;'];
    else
        output_file{end+1}=[slopeName '[' num2str(i-1) ']=' thisFormula ';'];
    end
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
        
        
        
        %convert "a^b" to "pow(a,b)" - find a(:),^,b(:) then replace with
        %{'pow' '(' a(:) ',' b(:) ')'}
        
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
                tok=tok(1:find(tok~='0',1,'last'));
                tok=[tok,'0'];
                
                if strcmpi(precision,'single') %add an f to indicate single precision
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