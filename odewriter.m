function odewriter(srcfilename, destinationFile)
% read ode write new ode (with options for sorting pars/vars/fixed/etc)

%parser already re-orders fixed quantities so they are always defined
%before use.

sortmethod='alpha';
% sortmethod='none'; %xpp is sensitive to order of equations...

fileExtension='ode';

% input checking
if nargin==0 || isempty(srcfilename)
    [name,path]=uigetfile('.ode','Select an ODE file');
    srcfilename=fullfile(path,name);
end

% Extract info from the ODE file
xppdata=parseODEfile(srcfilename);

%Build the destination filename
if ~exist('destinationFile','var')||isempty(destinationFile)
    %     destinationFile=xppdata.name;
    
    destinationFile=uiputfile('*.ode','New file name');
end

destinationFile=[destinationFile '.' fileExtension];
%end build destination filename

output_file=BuildOutputFile(xppdata,sortmethod);
lineCount=length(output_file);

%delete a previous version - silently destroys any previous version!
fullPath=[pwd filesep destinationFile];
if exist(fullPath,'file')==2
    delete(fullPath)
end

% write to file
fidw=fopen(fullPath,'w','n','UTF-8');

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

success = true;
return

end

function output_file=BuildOutputFile(xppdata,sortmethod)

% build the output as a cell array containing each line.
par=xppdata.par;
num=xppdata.num;
var=xppdata.var;
fixed=xppdata.fixed;
func=xppdata.func;
aux=xppdata.aux;
opt=xppdata.opt;


output_file={};
output_file{end+1}=['# ',xppdata.name];
output_file{end+1}='';

%constants
if xppdata.nNum>0
    output_file{end+1}='#Constants';
    num=sortNames(num,sortmethod);
    for i=1:xppdata.nNum
        output_file{end+1}=['num ' num(i).name '=' num2str(num(i).value)];
    end
end

%parameters
if xppdata.nPar>0
    output_file{end+1}='';
    output_file{end+1}='#Parameters';
    par=sortNames(par,sortmethod);
    for i=1:xppdata.nPar
        output_file{end+1}=['par ' par(i).name '=' num2str(par(i).value)];
    end
end

%Functions
if xppdata.nFunc>0
    output_file{end+1}='';
    output_file{end+1}='#Functions';
    func=sortNames(func,sortmethod);
    for i=1:xppdata.nFunc
        output_file{end+1}=[func(i).name '(' func(i).full_arglist ')=' func(i).formula];
    end
end

%Expressions
if xppdata.nFixed>0
    output_file{end+1}='';
    output_file{end+1}='#Fixed expressions';
    fixed=sortNames(fixed,sortmethod);
    for i=1:xppdata.nFixed
        output_file{end+1}=[fixed(i).name '=' fixed(i).formula];
    end
end


if xppdata.nVar>0
    var=sortNames(var,sortmethod);
    
    %Initial conditions
    output_file{end+1}='';
    output_file{end+1}='#Initial conditions';
    for i=1:xppdata.nVar
        output_file{end+1}=[var(i).name '(0)=' num2str(var(i).value)];
    end
    
    %State Variables
    output_file{end+1}='';
    output_file{end+1}='#Differential equations';
    for i=1:xppdata.nVar
        output_file{end+1}=[var(i).name '''=' var(i).formula];
    end
end

%User Functions
if xppdata.nAux>0
    output_file{end+1}='';
    output_file{end+1}='#User functions';
    aux=sortNames(aux,sortmethod);
    for i=1:xppdata.nAux
        output_file{end+1}=['aux ' aux(i).name '=' aux(i).formula];
    end
end

%options
if xppdata.nOpt>0
    output_file{end+1}='';
    output_file{end+1}='#Options';
    opt=sortNames(opt,sortmethod);
    for i=1:xppdata.nOpt
        output_file{end+1}=['@ ' opt(i).name '=' opt(i).value];
    end
end

output_file{end+1}='';
output_file{end+1}='done';

end

function mystruct=sortNames(mystruct,sortmethod)
switch sortmethod
    case 'alpha'
        [~,ixs]=sort({mystruct(:).name});
    otherwise
        ixs=1:numel(mystruct);
end
mystruct=mystruct(ixs);
end