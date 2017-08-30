function options=plotxppautset2(varargin)

% index of type of curve:
% 1-LP, 2-SNP, 3-HB, 4-TR?, 5-BR, 6-PD, 7-UZ/FP?

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('Seven entries per option correspond to: Saddle node, Saddle node of Periodics, Hopf, ?, Branch Point, Period Doubling, Fixed period \n');
    fprintf('           Color: [ {''k'',''k'',''r'',''b'',''k'',''r'',''b''} ]\n');
    fprintf('       LineStyle: [ {''-'',''-'',''none'',''none'',''-'',''none'',''none''} ]\n');
    fprintf('       Linewidth: [ {1, 1, 1, 1, 1, 1, 1} ]\n');
    fprintf('          Marker: [ {''none'',''none'',''none'',''none'',''none'',''none'',''none''} ]\n');
    fprintf('      Markersize: [ {2, 2, 3, 3, 2, 3, 3} ]\n');
    fprintf(' MarkerFaceColor: [ {''none'',''none'',''none'',''none'',''none'',''none'',''none''} ]\n');
    fprintf('      nSubSample: [1,1,1,1,1,1,1] \n')
    fprintf('\n');
    return;
end

Names = [
    'Color           '
    'LineStyle       '
    'Linewidth       '
    'Marker          '
    'Markersize      '
    'MarkerFaceColor '
    'nSubSample      '
    ];
m = size(Names,1);
names = lower(Names);



%defaults

%Color -- SHOULD BE VALID MATLAB Color
options.Color.SN = 'k';   %saddle node
options.Color.SNP = 'k';   %saddle node of periodics
options.Color.HB = 'r';   %hopf bifurcation
options.Color.uk1 = 'b';   %
options.Color.uk2 = 'b';   %
options.Color.uk3 = 'b';   %
options.Color.FP = 'b';   %fixed period

%LineStyle -- SHOULD BE VALID MATLAB LINESTYLES
options.LineStyle.SN = '-';
options.LineStyle.SNP = '-';
options.LineStyle.HB = '-';
options.LineStyle.uk1 = '-';
options.LineStyle.uk2 = '-';
options.LineStyle.uk3 = '-';
options.LineStyle.FP = '-';

%Linewidths
options.Linewidth.SN = 1;
options.Linewidth.SNP = 1;
options.Linewidth.HB = 1;
options.Linewidth.uk1 = 1;
options.Linewidth.uk2 = 1;
options.Linewidth.uk3 = 1;
options.Linewidth.FP = 1;

%Markersize
options.Marker.SN = 'none';
options.Marker.SNP = 'none';
options.Marker.HB = 'none';
options.Marker.uk1 = 'none';
options.Marker.uk2 = 'none';
options.Marker.uk3 = 'none';
options.Marker.FP = 'none';

%Markersize
options.Markersize.SN = 2;
options.Markersize.SNP = 2;
options.Markersize.HB = 2;
options.Markersize.uk1 = 2;
options.Markersize.uk2 = 2;
options.Markersize.uk3 = 2;
options.Markersize.FP = 2;

%MarkerFaceColor -- SHOULD BE VALID MATLAB Color
options.MarkerFaceColor.SN = 'none';
options.MarkerFaceColor.SNP = 'none';
options.MarkerFaceColor.HB = 'none';
options.MarkerFaceColor.uk1 = 'none';
options.MarkerFaceColor.uk2 = 'none';
options.MarkerFaceColor.uk3 = 'none';
options.MarkerFaceColor.FP = 'none';

options.nSubSample=[1,1,1,1,1,1,1];

%change any defaults with given inputs:

% A finite state machine to parse name-value pairs.
i=1;
if rem(nargin-i+1,2) ~= 0
    error(message('MATLAB:plotxppautset:ArgNameValueMismatch'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(message('MATLAB:plotxppautset:NoPropName', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(message('MATLAB:plotxppautset:InvalidPropName', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                matches = deblank(Names(j(1),:));
                for k = j(2:length(j))'
                    matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
                end
                error(message('MATLAB:plotxppautset:AmbiguousPropName',arg,matches));
            end
        end
        expectval = 1;                      % we expect a value next
        
    else
        if strcmpi(lowArg,'nSubSample')
            options.nSubSample=arg;
        else
            options.(deblank(Names(j,:))).SN = arg{1};
            options.(deblank(Names(j,:))).SNP = arg{2};
            options.(deblank(Names(j,:))).HB = arg{3};
            options.(deblank(Names(j,:))).uk1 = arg{4};
            options.(deblank(Names(j,:))).uk2 = arg{5};
            options.(deblank(Names(j,:))).uk3 = arg{6};
            options.(deblank(Names(j,:))).FP = arg{7};
        end
        expectval = 0;
        
    end
    i = i + 1;
end

if expectval
    error(message('MATLAB:plotxppautset:NoValueForProp', arg));
end

