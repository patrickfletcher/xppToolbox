function [lh,xnull,ynull]=plotnullclines(filename,xNullProp,yNullProp,options)
% Import and plot data from files saved from XPP's nullcline menu.
% Note: XPP outputs consecutive pairs of points to be connected as line
% segments.
%
% Line properties for each nullcline can be specified with optional input
% arguments using one of the following formats:
%
%    xNullProp = 'linespec'   ... e.g. 'r-o'
%    xNullProp = {'PropertyName','PropertyValue',...}   ???Doesn't work???
%    xNullProp.PropertyName=PropertyValue (available since R2014b)
%
% Default format: no markers, x-nullcline orange, y-nullcline green.
%
% See line properties for full description of options available.
%
% (c) Patrick Fletcher 2017


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
else
    %update the full option structure with an "optionset" routine
end

%decode LineSpec type option input ->   %    xNullProp = 'LineSpec';
% if ischar(xNullProp)
%     xNullProp=decodeLineSpec(xNullProp);
% end
% if ischar(yNullProp)
%     yNullProp=decodeLineSpec(yNullProp)
% end

[xnull, ynull]=readNullclineFile(filename);

% xnull=xnull(1:xOpts.subsample:end,:);
% ynull=ynull(1:yOpts.subsample:end,:);

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

% %line segment method - only looks nice with solid line format
% hold on
% lh.x=[];
% for i=1:2:nX
%     lh.x(end+1)=plot([xnull(i,1) xnull(i+1,1)],[xnull(i,2) xnull(i+1,2)],xNullProp);
% %     lh.x(end+1)=line([xnull(i,1) xnull(i+1,1)],[xnull(i,2) xnull(i+1,2)],xNullProp);
% end
%
% lh.y=[];
% for i=1:2:nY
%     lh.y(end+1)=plot([ynull(i,1) ynull(i+1,1)],[ynull(i,2) ynull(i+1,2)],yNullProp);
% %     lh.y(end+1)=line([ynull(i,1) ynull(i+1,1)],[ynull(i,2) ynull(i+1,2)],yNullProp);
% end
% hold off

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


function props=decodeLineSpec(str)
% Decode linespec string str, convert to property/value struct
%
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot
%          c     cyan          +     plus               --    dashed
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%          w     white         v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram

color={'b','g','r','c','m','y','k','w'};
marker={'.','o','x','+','*','s','d','v','^','<','>','p','h'};
linestyle={'-',':','-.','--'};

foundColor=false;
foundMarker=false;
foundLinestyle=false;

i=1;
while i<=length(str)
    
    isColor=strcmpi(color,str(i));
    isMarker=strcmpi(marker,str(i));
    isLinestyle=strcmpi(linestyle,str(i));
    
    if any(isColor) && ~foundColor
        props.Color=color{isColor};
        foundColor=true;
        
    elseif any(isMarker) && ~foundMarker
        props.Marker=marker{isMarker};
        foundMarker=true;
        
    elseif any(isLinestyle) && ~foundLinestyle
        props.LineStyle=linestyle{isLinestyle}; %'-', ':'
        
        if strcmpi(str(i),'-') && i<length(str) %must check if '-.' or '--'
            next=str(i+1);
            if strcmpi(next,'-')
                props.LineStyle=[props.LineStyle '-'];
                i=i+1; %skip ahead
            elseif strcmpi(next,'.')
                props.LineStyle=[props.LineStyle '.'];
                i=i+1;
            end
        end
        foundLinestyle=true;
    else
        error('Error in color/linetype argument')
    end
    
    i=i+1;
end

%default output to allow consistency with plot
if ~foundLinestyle
    if foundMarker
        props.LineStyle='none'; %if marker specified without linestyle, assume linestyle is none
    else
        props.LineStyle='-'; %no marker, no linestyle => use '-' and optionally specified color
    end
end

%if desired, could return cell containing name/value pair instead...
% props=struct2nameValue(props);

end
