function plotxppaut1(filename,options)
% (c) Patrick Fletcher 2015
% Based on PlotXppaut by Mohammad S. Imtiaz
%
% A function to plot one-parameter bifurcation diagrams from XPPAUT.
% Provides ability to set options to customize all aspects of the plot.
% downsampling of points in the dataset for plotting and externally set
% linespecs

% 07/08/2016: Support XPP8.0's 6 column format. The extra column is to
% specify twoparameter bifurcation type. Two parameter bifurcation diagrams
% will require a different way to map linespecs to curves?

%file format - columns: V>=8.0 {par, ymin, ymax, type, branch, 2-par type};
% V<8.0 {par, ymin, ymax, type, branch}

nSubSample=options.nSubSample;

if ~exist('options','var')
    options=plotxppautset1();
end

%Colors -- SHOULD BE VALID MATLAB COLORS
C=options.Colors;

%LineStyle -- SHOULD BE VALID MATLAB LINESTYLES
Lt=options.LineStyle;

%Linewidths
Lw=options.Linewidths;

%Marker
M=options.Marker;

%Markersize
Ms=options.Markersize;

%MarkerFaceColor -- SHOULD BE VALID MATLAB COLORS
Mfc=options.MarkerFaceColor;


fid = fopen(filename,'rt');

if(fid~=-1)
    
    %determine if old or new format
    %src: https://www.mathworks.com/matlabcentral/answers/54236-find-number-of-columns-on-text-file
    
    delimiter = ' '; %or whatever
    tLines = fgets(fid);
    numCols = numel(strfind(tLines,delimiter)) + 1;
    
    st = fscanf(fid,'%f',[numCols,inf]); %in Xpp version 8 - col 6 stores type of Two-par curves
    fclose(fid);
    
    
    temp=st';
    
    temp1=zeros(size(temp,1),8)*NaN;
    
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
    
    temp1(type1ix,[1 5]) = temp(type1ix,[2 3]);
    temp1(type2ix,[2 6]) = temp(type2ix,[2 3]);
    temp1(type3ix,[3 7]) = temp(type3ix,[2 3]);
    temp1(type4ix,[4 8]) = temp(type4ix,[2 3]);
    
    
    %above does the same without loops.
    %   for n=1:size(temp1,1),
    %     if(temp(n,4)==1),
    %       temp1(n,[1 5]) = temp(n,[2 3]);
    %     end;
    %     if(temp(n,4)==2),
    %       temp1(n,[2 6]) = temp(n,[2 3]);
    %     end;
    %     if(temp(n,4)==3),
    %       temp1(n,[3 7]) = temp(n,[2 3]);
    %     end;
    %     if(temp(n,4)==4),
    %       temp1(n,[4 8]) = temp(n,[2 3]);
    %     end;
    %   end;
    
    
    % $$$ figure;
    %	-----------------------------------------
    h_SS = plot(temp(:,1),temp1(:,[1 5]));
    set(h_SS,'color',C.Ss)
    set(h_SS,'linestyle',Lt.Ss)
    set(h_SS,'linewidth',Lw.Ss)
    set(h_SS,'Marker',M.Ss)
    set(h_SS,'Markersize',Ms.Ss)
    set(h_SS,'MarkerFacecolor',Mfc.Ss)
    hold on;
    
    %	-----------------------------------------
    h_US = plot(temp(:,1),temp1(:,[2 6]));
    set(h_US,'color',C.Us)
    set(h_US,'linestyle',Lt.Us)
    set(h_US,'linewidth',Lw.Us)
    set(h_US,'Marker',M.Us)
    set(h_US,'Markersize',Ms.Us)
    set(h_US,'MarkerFacecolor',Mfc.Us)
    
    %	-----------------------------------------
    
    
    h_SP = plot(temp(:,1),temp1(:,[3 7]));
    set(h_SP,'color',C.Sp)
    set(h_SP,'linestyle',Lt.Sp)
    set(h_SP,'linewidth',Lw.Sp)
    set(h_SP,'Marker',M.Sp)
    set(h_SP,'Markersize',Ms.Sp)
    set(h_SP,'MarkerFacecolor',Mfc.Sp)
    %	-----------------------------------------
    
    
    h_UP = plot(temp(:,1),temp1(:,[4 8]));
    set(h_UP,'color',C.Up)
    set(h_UP,'linestyle',Lt.Up)
    set(h_UP,'linewidth',Lw.Up)
    set(h_UP,'Marker',M.Up)
    set(h_UP,'Markersize',Ms.Up)
    set(h_UP,'MarkerFacecolor',Mfc.Up)
    %	-----------------------------------------
    
    
    %   grid on;
    %   axis tight
    %   xlabel('Bifurcation parameter')
    %   ylabel('Variable')
    %   title(file_name);
else
    error(['Could not open file: ' filename ])
end
