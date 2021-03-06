function success=ode2vcml(sourceFilename, destinationFilename)
%ODE2VCML Convert an ODE file for XPP to a VCML file for Virtual Cell
% The generated .vcml file can be imported as a Math Model.
%
%USAGE:
% ODE2VCML() prompts the user for .ode file, writes to a file of the same
% name but with .vcml extension
%
% ODE2VCML(sourceFilename) reads specified .ode file
%
% ODE2VCML(sourceFilename,destinationFilename) reads specified .ode file,
% and writes to specified destination .vcml file
%
%NOTES:
% A subset of XPP syntax is supported:
%
%       Parameters:
%           p* name=value, name=value, ...
%
%       Constants:
%           n* name=value, name=value, ...
%
%       Ordinary differential equations:
%           name'=formula
%           dnamedt=formula
%
%       Initial conditions:
%           name(0)=value
%           init name=value, name=value, ...
%
%       Fixed expressions (temporary quantities)
%           name=formula
%           !name = formula
%
%       Auxiliary variables:
%           aux name=formula
%
%       XPP options:
%           @ name=value, name=value, ...
%
%       Comments:
%           # comment_text
%           % comment_text
%           " comment_text
%

% Copyright 2016 Patrick Allen Fletcher
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%get an ode file if needed
if nargin==0 || isempty(sourceFilename)
    [name,path]=uigetfile('.ode','Select an ODE file');
    sourceFilename=fullfile(path,name);
end

% Extract all info from the ODE file
xppdata=parseODEfile(sourceFilename);

%Destination filename
if ~exist('destinationFilename','var')||isempty(destinationFilename)
    destinationFilename=xppdata.name;
end
if ~contains(destinationFilename,'.vcml')
    destinationFilename=[destinationFilename '.vcml'];
end

%Destination file exists?
fullPath=[pwd filesep destinationFilename];
if exist(fullPath,'file')==2
    choice = questdlg('Destination file exists. Overwrite?', ...
        'Yes','No');
    % Handle response
    switch choice
        case 'Yes'
            delete(fullPath)
        case 'No'
            [filename, pathname] = uiputfile('*.vcml','Save as');
            fullPath=[pathname,filename];
        case 'Cancel'
            error('User canceled - aborting...')
    end
end

%build the output file contents
output_file=BuildOutputFile(xppdata);

% write to file
fidw=fopen(fullPath,'w','n','UTF-8');

lineCount=length(output_file);
if fidw > -1
    for line=1:lineCount
        fprintf(fidw, '%s\n', output_file{line});
    end
    status=fclose(fidw);
else
    fprintf(' Problem opening file %s for writing. Cannot continue\n',fullPath);
    return
end

success=true;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output_file=BuildOutputFile(xppdata)


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