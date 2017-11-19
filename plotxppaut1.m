function h=plotxppaut1(filename,options)
% plotxppaut2 Plot one-parameter bifurcation diagrams generated in XPPAUT.
% Generate the data file data file using XPPAUT's Auto/File/write points
% option.
%
%   plotxppaut1(filename) - plot contents of filename
%   plotxppaut1(filename,options) - plot with line properties defined in
%   the options struct, which can be generated and modified using plotxppautset2
%
% See also  plotxppautset1 plotxppaut2 plotxppautset2 plotnullclines

% (c) Patrick Fletcher 2015
% Inspired by PlotXppaut by Mohammad S. Imtiaz
% 07/08/2016: Support XPP8.0's 6 column format. The extra column is to
% specify two-parameter bifurcation type. 
%
% columns: V>=8.0 {par, ymin, ymax, type, branch, 2-par type};
% V<8.0 {par, ymin, ymax, type, branch}

%TODO: use struct of arrays for lineprops -> index easily via type
%TODO: plot by branch not type?
%TODO: handle change in type: no gaps in lines..

if ~exist('options','var')
    options=plotxppautset1();
end
nSubSample=options.nSubSample;

%Colors -- SHOULD BE VALID MATLAB COLORS
C=options.Color;

%LineStyle -- SHOULD BE VALID MATLAB LINESTYLES
Lt=options.LineStyle;

%Linewidths
Lw=options.Linewidth;

%Marker
M=options.Marker;

%Markersize
Ms=options.Markersize;

%MarkerFaceColor -- SHOULD BE VALID MATLAB COLORS
Mfc=options.MarkerFaceColor;

fid = fopen(filename,'rt');

h=[];

if(fid~=-1)
    
    %determine if old or new format
    %src: https://www.mathworks.com/matlabcentral/answers/54236-find-number-of-columns-on-text-file
    
    delimiter = ' '; %or whatever
    tLines = fgets(fid);
    numCols = numel(strfind(tLines,delimiter)) + 1;
    
    st = fscanf(fid,'%f',[numCols,inf]); %in Xpp version 8 - col 6 stores type of Two-par curves
    fclose(fid);
    
    
    temp=st';
    
    x=nan(size(temp,1),4);
    y=nan(size(temp,1),4);
    ymin=nan(size(temp,1),2);
    ymax=nan(size(temp,1),2);
    
    type1=temp(:,4)==1;
    type2=temp(:,4)==2;
    type3=temp(:,4)==3;
    type4=temp(:,4)==4;
    
    type1ix=find(type1);
    type2ix=find(type2);
    type3ix=find(type3);
    type4ix=find(type4);
    
    %subsample
    type1ix=type1ix(1:nSubSample(1):end);
    type2ix=type2ix(1:nSubSample(2):end);
    type3ix=type3ix(1:nSubSample(3):end);
    type4ix=type4ix(1:nSubSample(4):end);
    
    x(type1ix,1)=temp(type1ix,1);
    x(type2ix,2)=temp(type2ix,1);
    x(type3ix,3)=temp(type3ix,1);
    x(type4ix,4)=temp(type4ix,1);
    
    y(type1ix,1)=temp(type1ix,2);
    y(type2ix,2)=temp(type2ix,2);
    ymin(type3ix,1)=temp(type3ix,2);
    ymax(type3ix,1)=temp(type3ix,3);
    ymin(type4ix,2)=temp(type4ix,2);
    ymax(type4ix,2)=temp(type4ix,3);
    
%     temp1(type1ix,[1 5]) = temp(type1ix,[2 3]);
%     temp1(type2ix,[2 6]) = temp(type2ix,[2 3]);
%     temp1(type3ix,[3 7]) = temp(type3ix,[2 3]);
%     temp1(type4ix,[4 8]) = temp(type4ix,[2 3]);
    
    
    %%% Now plot
    %	-----------------------------------------
%     h.SS = plot(temp(:,1),temp1(:,[1 5]));
    h.SS = plot(x(:,1),y(:,1));
    set(h.SS,'color',C.Ss)
    set(h.SS,'linestyle',Lt.Ss)
    set(h.SS,'linewidth',Lw.Ss)
    set(h.SS,'Marker',M.Ss)
    set(h.SS,'Markersize',Ms.Ss)
    set(h.SS,'MarkerFacecolor',Mfc.Ss)
    hold on;
    
    %	-----------------------------------------
%     h.US = plot(temp(:,1),temp1(:,[2 6]));
    h.US = plot(x(:,2),y(:,2));
    set(h.US,'color',C.Us)
    set(h.US,'linestyle',Lt.Us)
    set(h.US,'linewidth',Lw.Us)
    set(h.US,'Marker',M.Us)
    set(h.US,'Markersize',Ms.Us)
    set(h.US,'MarkerFacecolor',Mfc.Us)
    
    %	-----------------------------------------
    
    
%     h.SP = plot(temp(:,1),temp1(:,[3 7]));
    h.SP = plot(x(:,3),[ymin(:,1), ymax(:,1)]);
    set(h.SP,'color',C.Sp)
    set(h.SP,'linestyle',Lt.Sp)
    set(h.SP,'linewidth',Lw.Sp)
    set(h.SP,'Marker',M.Sp)
    set(h.SP,'Markersize',Ms.Sp)
    set(h.SP,'MarkerFacecolor',Mfc.Sp)
    %	-----------------------------------------
    
    
%     h.UP = plot(temp(:,1),temp1(:,[4 8]));
    h.UP = plot(x(:,4),[ymin(:,2), ymax(:,2)]);
    set(h.UP,'color',C.Up)
    set(h.UP,'linestyle',Lt.Up)
    set(h.UP,'linewidth',Lw.Up)
    set(h.UP,'Marker',M.Up)
    set(h.UP,'Markersize',Ms.Up)
    set(h.UP,'MarkerFacecolor',Mfc.Up)
    %	-----------------------------------------
    
    
else
    error(['Could not open file: ' filename ])
end
