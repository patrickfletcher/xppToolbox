function [lh,xnull,ynull]=plotnullclines(filename,xNullProp,yNullProp,options)
% PLOTNULLCLINES Import and plot data from files saved from XPP's nullcline menu.
% Note: XPP outputs consecutive pairs of points to be connected as line
% segments.
%
%   plotnullclines(filename) - plot contents of filename using default colors
%
%   plotnullclines(filename,xNullProp,yNullProp) - plot using line
%   properties specified in xNullProp and/or yNullProp. To use default,
%   specify as [].
%
%   plotnullclines(filename,xNullProp,yNullProp,options) - specify how to
%   sort the x and y nullcline data in a struct:
%       option.sortXByX={true}/false - sort x-nullcline by x values
%       option.sortYByX={true}/false - sort y-nullcline by x values
%       options.XkeepEvery = {1} - for subsampling
%       options.YkeepEvery = {1} - for subsampling
%
%   [lh,xnull,ynull]=plotnullclines(filename,...) - return line handles and
%   nullcline data
%
% Line properties for each nullcline can be specified with optional input
% arguments using one of the following formats:
%
%    xNullProp = 'linespec'   ... e.g. 'r-o'
%    xNullProp.PropertyName=PropertyValue (available since R2014b)
%
% Default format: no markers, x-nullcline orange, y-nullcline green. See line properties for full description of options available.
%
%    options
%
% See also  plotxppautset1 plotxppaut2 plotxppautset2 plotnullclines

% (c) Patrick Fletcher 2017

% TODO: best way, struct or spec?   xNullProp = {'PropertyName','PropertyValue',...}   ???Doesn't work???

%defaults that look like xpp.
if ~exist('xNullProp','var') || isempty(xNullProp)
    xNullProp.Color=[0.75,0.5,0];
end
if ~exist('yNullProp','var') || isempty(yNullProp)
    yNullProp.Color=[0,0.75,0];
end
if ~exist('options','var') || isempty(options)
    options.sortXByX=true;
    options.sortYByX=true;
    options.XkeepEvery=1;
    options.YkeepEvery=1;
else
    %update the full option structure with an "optionset" routine
end

[xnull, ynull]=readNullclineFile(filename);

xnull=xnull(1:options.XkeepEvery:end,:);
ynull=ynull(1:options.YkeepEvery:end,:);

%mimic behavior of hold on/off from "plot"
gcf;
if ~ishold
    cla
end

%sort and plot method (must be a proper function of one of x or y...)
if options.sortXByX
    [~,ix]=sort(xnull(:,1)); %by x-values
else
    [~,ix]=sort(xnull(:,2)); %by y-values
end
xnull=xnull(ix,:); 

if options.sortYByX
    [~,iy]=sort(ynull(:,1)); %by x-values
else
    [~,iy]=sort(ynull(:,2)); %by y-values
end
ynull=ynull(iy,:);

hold on
lh.x=plot(xnull(:,1),xnull(:,2),xNullProp);
lh.y=plot(ynull(:,1),ynull(:,2),yNullProp);
hold off

end

function [xnull, ynull]=readNullclineFile(filename)
% reads the dat file. Interprets 'x y 1' and 'x y 2' as the (x,y)
% coordinates for the X and Y nullclines. Ignores all other lines in the
% file.

xnull=[];
ynull=[];

fid=fopen(filename);
fline = fgetl(fid);
while ischar(fline)
    
    if ~isempty(fline)
        C=textscan(fline,'%f');
    else
        C{1}=[];
    end
    
    if ~isempty(C{1})
        if C{1}(3)==1 %xnull
            xnull(end+1,1)=C{1}(1);
            xnull(end,2)=C{1}(2);
        elseif C{1}(3)==2 %ynull
            ynull(end+1,1)=C{1}(1);
            ynull(end,2)=C{1}(2);
        end
    end
    
    fline = fgetl(fid);
end
fclose(fid);
end