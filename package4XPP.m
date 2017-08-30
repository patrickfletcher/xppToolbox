function xppStruct=package4XPP(pNames,pVals,yNames,y0Vals)
% PACKAGE4XPP Package parameters and/or initial conditions into a format
% for use with ChangeXPPodeFile.
%
%   xppStruct=PACKAGE4XPP(pNames,pVals) - parameters only
%   xppStruct=PACKAGE4XPP([],[],yNames,y0Vals) - initial conditions only
%   xppStruct=PACKAGE4XPP(pNames,pVals,yNames,y0Vals)
%
% See also ChangeXPPodeFile, ChangeXPPsetFile

%TODO: rename?
%TODO: update ChangeXPPodeFile, ChangeXPPsetFile (link to original)

if nargin<2||nargin==3||nargin>4
    error('incorrect number of arguments')
end
if nargin>1
    nPar=length(pVals);
    nY0=0;
end
if nargin>3
    nY0=length(y0Vals);
end

nElems=nPar+nY0;

xppStruct=struct('type',cell(nElems,1),'name',cell(nElems,1),'val',cell(nElems,1));

for i=1:nPar
    xppStruct(i).type='PAR';
    xppStruct(i).name=pNames{i};
    xppStruct(i).val=pVals(i);
end

for i=1:nY0
    xppStruct(nPar+i).type='IC';
    xppStruct(nPar+i).name=yNames{i};
    xppStruct(nPar+i).val=y0Vals(i);
end