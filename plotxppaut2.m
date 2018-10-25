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

%TODO: use struct of arrays for lineprops -> index easily via type
%TODO: handle multiple parts of same branch eg. ds<0 then ds>0??? large
%jumps in dx,dy?
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

% h=[];

if(fid~=-1)
    
    %determine if old or new format
    %src: https://www.mathworks.com/matlabcentral/answers/54236-find-number-of-columns-on-text-file
    
    delimiter = ' '; %or whatever
    tLines = fgets(fid);
    numCols = numel(strfind(tLines,delimiter)) + 1;
    
    st = fscanf(fid,'%f',[numCols,inf]); %in Xpp version 8 - col 6 stores type of Two-par curves
    fclose(fid);
    
    
    temp=st';
    branch=[temp(:,5),temp(:,6)];
    [ubr, ubrix]=unique(branch,'rows','stable');
    
    twopartype=temp(:,6);
    
    %%% plot each branch, line properties set according to type
    for i=1:size(ubr,1)
        thistype=twopartype(ubrix(i));
        thisbr=branch(:,1)==ubr(i,1)&branch(:,2)==ubr(i,2);
        thisbrix=find(thisbr);
        
        xx=temp(thisbrix(1:nSubSample(thistype):end),1);
        yy=temp(thisbrix(1:nSubSample(thistype):end),2);
        
        h(i)=line(xx,yy);
        
        setLineProps(h(i),thistype)
        
    end
    
    
else
    error(['Could not open file: ' filename ])
end

function setLineProps(line_handle,thisType)
    switch thisType
        case 1
            set(line_handle,'color',C.SN)
            set(line_handle,'linestyle',Lt.SN)
            set(line_handle,'linewidth',Lw.SN)
            set(line_handle,'Marker',M.SN)
            set(line_handle,'Markersize',Ms.SN)
            set(line_handle,'MarkerFacecolor',Mfc.SN)
        case 2
            set(line_handle,'color',C.SNP)
            set(line_handle,'linestyle',Lt.SNP)
            set(line_handle,'linewidth',Lw.SNP)
            set(line_handle,'Marker',M.SNP)
            set(line_handle,'Markersize',Ms.SNP)
            set(line_handle,'MarkerFacecolor',Mfc.SNP)
        case 3
            set(line_handle,'color',C.HB)
            set(line_handle,'linestyle',Lt.HB)
            set(line_handle,'linewidth',Lw.HB)
            set(line_handle,'Marker',M.HB)
            set(line_handle,'Markersize',Ms.HB)
            set(line_handle,'MarkerFacecolor',Mfc.HB)
        case 4
            set(line_handle,'color',C.uk1)
            set(line_handle,'linestyle',Lt.uk1)
            set(line_handle,'linewidth',Lw.uk1)
            set(line_handle,'Marker',M.uk1)
            set(line_handle,'Markersize',Ms.uk1)
            set(line_handle,'MarkerFacecolor',Mfc.uk1)
        case 5
            set(line_handle,'color',C.uk2)
            set(line_handle,'linestyle',Lt.uk2)
            set(line_handle,'linewidth',Lw.uk2)
            set(line_handle,'Marker',M.uk2)
            set(line_handle,'Markersize',Ms.uk2)
            set(line_handle,'MarkerFacecolor',Mfc.uk2)   
        case 6
            set(line_handle,'color',C.uk3)
            set(line_handle,'linestyle',Lt.uk3)
            set(line_handle,'linewidth',Lw.uk3)
            set(line_handle,'Marker',M.uk3)
            set(line_handle,'Markersize',Ms.uk3)
            set(line_handle,'MarkerFacecolor',Mfc.uk3)
        case 7
            set(line_handle,'color',C.FP)
            set(line_handle,'linestyle',Lt.FP)
            set(line_handle,'linewidth',Lw.FP)
            set(line_handle,'Marker',M.FP)
            set(line_handle,'Markersize',Ms.FP)
            set(line_handle,'MarkerFacecolor',Mfc.FP)
    end
end

end