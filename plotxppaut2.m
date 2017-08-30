function h=plotxppaut2(filename,options)
% plotxppaut2 Plot two-parameter bifurcation diagrams generated in XPPAUT.
% Generate the data file data file using XPPAUT's Auto/File/write points
% option.
%
%   plotxppaut2(filename) - plot contents of filename
%   plotxppaut2(filename,options) - plot with line properties defined in
%   the options struct, which can be generated and modified using plotxppautset2
%
% See also  plotxppautset2 plotxppaut1 plotxppautset1 plotnullclines

% (c) Patrick Fletcher 2015-now
% Inspired by PlotXppaut by Mohammad S. Imtiaz
%
% 07/08/2016: Support XPP8.0's 6 column format. The extra column is to
% specify two-parameter bifurcation type.
%
% columns: V>=8.0 {par, ymin, ymax, type, branch, 2-par type};
% V<8.0 {par, ymin, ymax, type, branch}

%TODO: confirm point types
%TODO: MORE INTUITIVE options struct, for each type of branch, set linespec - e.g. 'SNP','-ro'

if ~exist('options','var')
    options=plotxppautset2();
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
    
    type1=temp(:,6)==1;
    type2=temp(:,6)==2;
    type3=temp(:,6)==3;
    type4=temp(:,6)==4;
    type5=temp(:,6)==5;
    type6=temp(:,6)==6;
    type7=temp(:,6)==7;
    
    type1ix=find(type1);
    type2ix=find(type2);
    type3ix=find(type3);
    type4ix=find(type4);
    type5ix=find(type5);
    type6ix=find(type6);
    type7ix=find(type7);
    
    %subsample
    type1ix=type1ix(1:nSubSample(1):end);
    type2ix=type2ix(1:nSubSample(2):end);
    type3ix=type3ix(1:nSubSample(3):end);
    type4ix=type4ix(1:nSubSample(4):end);
    type5ix=type5ix(1:nSubSample(5):end);
    type6ix=type6ix(1:nSubSample(6):end);
    type7ix=type7ix(1:nSubSample(7):end);
    
    x(type1ix,1)=temp(type1ix,1);
    x(type2ix,2)=temp(type2ix,1);
    x(type3ix,3)=temp(type3ix,1);
    x(type4ix,4)=temp(type4ix,1);
    x(type5ix,5)=temp(type5ix,1);
    x(type6ix,6)=temp(type6ix,1);
    x(type7ix,7)=temp(type7ix,1);
    
    y(type1ix,1)=temp(type1ix,2);
    y(type2ix,2)=temp(type2ix,2);
    y(type3ix,3)=temp(type3ix,2);
    y(type4ix,4)=temp(type4ix,2);
    y(type5ix,5)=temp(type5ix,2);
    y(type6ix,6)=temp(type6ix,2);
    y(type7ix,7)=temp(type7ix,2);
    
    
    
    %%% Now plot
    %	-----------------------------------------
    h.one = line(x(:,1),y(:,1));
    set(h.one,'color',C.SN)
    set(h.one,'linestyle',Lt.SN)
    set(h.one,'linewidth',Lw.SN)
    set(h.one,'Marker',M.SN)
    set(h.one,'Markersize',Ms.SN)
    set(h.one,'MarkerFacecolor',Mfc.SN)
    
    %	-----------------------------------------
    h.two = line(x(:,2),y(:,2));
    set(h.two,'color',C.SNP)
    set(h.two,'linestyle',Lt.SNP)
    set(h.two,'linewidth',Lw.SNP)
    set(h.two,'Marker',M.SNP)
    set(h.two,'Markersize',Ms.SNP)
    set(h.two,'MarkerFacecolor',Mfc.SNP)
    
    %	-----------------------------------------
    h_three = line(x(:,3),y(:,3));
    set(h_three,'color',C.HB)
    set(h_three,'linestyle',Lt.HB)
    set(h_three,'linewidth',Lw.HB)
    set(h_three,'Marker',M.HB)
    set(h_three,'Markersize',Ms.HB)
    set(h_three,'MarkerFacecolor',Mfc.HB)
    
    %	-----------------------------------------
    h.four = line(x(:,4),y(:,4));
    set(h.four,'color',C.uk1)
    set(h.four,'linestyle',Lt.uk1)
    set(h.four,'linewidth',Lw.uk1)
    set(h.four,'Marker',M.uk1)
    set(h.four,'Markersize',Ms.uk1)
    set(h.four,'MarkerFacecolor',Mfc.uk1)
    
    %	-----------------------------------------
    h.five = line(x(:,5),y(:,5));
    set(h.five,'color',C.uk2)
    set(h.five,'linestyle',Lt.uk2)
    set(h.five,'linewidth',Lw.uk2)
    set(h.five,'Marker',M.uk2)
    set(h.five,'Markersize',Ms.uk2)
    set(h.five,'MarkerFacecolor',Mfc.uk2)
    
    %	-----------------------------------------
    h.six = line(x(:,6),y(:,6));
    set(h.six,'color',C.uk3)
    set(h.six,'linestyle',Lt.uk3)
    set(h.six,'linewidth',Lw.uk3)
    set(h.six,'Marker',M.uk3)
    set(h.six,'Markersize',Ms.uk3)
    set(h.six,'MarkerFacecolor',Mfc.uk3)
    
    %	-----------------------------------------
    h.seven = line(x(:,7),y(:,7));
    set(h.seven,'color',C.FP)
    set(h.seven,'linestyle',Lt.FP)
    set(h.seven,'linewidth',Lw.FP)
    set(h.seven,'Marker',M.FP)
    set(h.seven,'Markersize',Ms.FP)
    set(h.seven,'MarkerFacecolor',Mfc.FP)
    
    %	-----------------------------------------
    
else
    error(['Could not open file: ' filename ])
end
