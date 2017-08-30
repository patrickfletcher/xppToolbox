function plotxppaut2(filename,options)
% Plot two parameter bifurcation diagrams using the new XPP8.0 write points
% format. Column 6 contains an indicator of the type of bifurcation point,
% so we can use that to assign linespecs.

% Need to discover the ID code:
% torus, bifurcations
%
% MORE INTUITIVE: for each type of branch, set linespec - e.g. 'SNP','-ro'


nSubSample=options.nSubSample;

if ~exist('options','var')
    options=plotxppautset2();
end

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

if(fid~=-1)
    
    %determine if old or new format
    %src: https://www.mathworks.com/matlabcentral/answers/54236-find-number-of-columns-on-text-file
    
    delimiter = ' '; %or whatever
    tLines = fgets(fid);
    numCols = numel(strfind(tLines,delimiter)) + 1;
    
    st = fscanf(fid,'%f',[numCols,inf]); %in Xpp version 8 - col 6 stores type of Two-par curves
    fclose(fid);
    
    
    temp=st';
    
    temp1=zeros(size(temp,1),7)*NaN;
    
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
    %   type1ix=type1ix(1:nSubSample(1):end);
    %   type2ix=type2ix(1:nSubSample(2):end);
    %   type3ix=type3ix(1:nSubSample(3):end);
    %   type4ix=type4ix(1:nSubSample(4):end);
    %   type5ix=type5ix(1:nSubSample(5):end);
    %   type6ix=type6ix(1:nSubSample(6):end);
    %   type7ix=type7ix(1:nSubSample(7):end);
    
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
    
    
    
    % $$$ figure;
    %	-----------------------------------------
    h_one = line(x(:,1),y(:,1));
    set(h_one,'color',C.SN)
    set(h_one,'linestyle',Lt.SN)
    set(h_one,'linewidth',Lw.SN)
    set(h_one,'Marker',M.SN)
    set(h_one,'Markersize',Ms.SN)
    set(h_one,'MarkerFacecolor',Mfc.SN)
    
    %	-----------------------------------------
    h_two = line(x(:,2),y(:,2));
    set(h_two,'color',C.SNP)
    set(h_two,'linestyle',Lt.SNP)
    set(h_two,'linewidth',Lw.SNP)
    set(h_two,'Marker',M.SNP)
    set(h_two,'Markersize',Ms.SNP)
    set(h_two,'MarkerFacecolor',Mfc.SNP)
    
    %	-----------------------------------------
    h_three = line(x(:,3),y(:,3));
    set(h_three,'color',C.HB)
    set(h_three,'linestyle',Lt.HB)
    set(h_three,'linewidth',Lw.HB)
    set(h_three,'Marker',M.HB)
    set(h_three,'Markersize',Ms.HB)
    set(h_three,'MarkerFacecolor',Mfc.HB)
    
    %	-----------------------------------------
    h_four = line(x(:,4),y(:,4));
    set(h_four,'color',C.uk1)
    set(h_four,'linestyle',Lt.uk1)
    set(h_four,'linewidth',Lw.uk1)
    set(h_four,'Marker',M.uk1)
    set(h_four,'Markersize',Ms.uk1)
    set(h_four,'MarkerFacecolor',Mfc.uk1)
    
    %	-----------------------------------------
    h_five = line(x(:,5),y(:,5));
    set(h_five,'color',C.uk2)
    set(h_five,'linestyle',Lt.uk2)
    set(h_five,'linewidth',Lw.uk2)
    set(h_five,'Marker',M.uk2)
    set(h_five,'Markersize',Ms.uk2)
    set(h_five,'MarkerFacecolor',Mfc.uk2)
    
    %	-----------------------------------------
    h_six = line(x(:,6),y(:,6));
    set(h_six,'color',C.uk3)
    set(h_six,'linestyle',Lt.uk3)
    set(h_six,'linewidth',Lw.uk3)
    set(h_six,'Marker',M.uk3)
    set(h_six,'Markersize',Ms.uk3)
    set(h_six,'MarkerFacecolor',Mfc.uk3)
    
    %	-----------------------------------------
    h_seven = line(x(:,7),y(:,7));
    set(h_seven,'color',C.FP)
    set(h_seven,'linestyle',Lt.FP)
    set(h_seven,'linewidth',Lw.FP)
    set(h_seven,'Marker',M.FP)
    set(h_seven,'Markersize',Ms.FP)
    set(h_seven,'MarkerFacecolor',Mfc.FP)
    
    %	-----------------------------------------
    
    %   grid on;
    %   axis tight
    %   xlabel('Bifurcation parameter')
    %   ylabel('Variable')
    %   title(file_name);
else
    error(['Could not open file: ' filename ])
end
