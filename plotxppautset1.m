function options=plotxppautset1(varargin)
% PLOTXPPAUTSET1 Create or update an option structure to specify the line properties of each
% type of solution branch in a 1-parameter bifurcation diagram to be
% plotted by plotxppaut1. Defaults are set to mimic XPPAUT appearance.
%
%    optionStruct=PLOTXPPAUTSET1('Name','Value')
%
% Where 'Name' and 'Value' are line properties; 
%
% See also  plotxppaut1 plotxppaut2 plotxppautset2 plotnullclines


% index of type of curve:
% 1-Stable SS, 2-Unstable SS, 3-Stable LC, 4-Unstable LC


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('Create or update an options structure with \n');
    fprintf('the following fields (defaults shown):\n');
    fprintf('\n');
    fprintf('          Color: {''k'',''k'',''r'',''b''}\n');
    fprintf('       LineStyle: {''-'',''--'',''none'',''none''}\n');
    fprintf('      Linewidth: {1.5, 1, 1, 1}\n');
    fprintf('          Marker: {''none'',''none'',''.'',''o''}\n');
    fprintf('      Markersize: {2, 2, 3, 3}\n');
    fprintf(' MarkerFaceColor: {''none'',''none'',''r'',''none''}\n');
    fprintf('      nSubSample: [1,1,1,1]\n')
    fprintf('\n');
    fprintf('The four entries per option correspond to:\n');
    fprintf('{Stable SS, Unstable SS, Stable LC, Unstable LC}\n');
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

%Colors -- SHOULD BE VALID MATLAB COLORS
options.Color.Ss = 'r';   %STABLE STEADY STATE
options.Color.Us = 'k';   %UNSTABLE STEADY STATE
options.Color.Sp = 'g';     %STABLE PERIODIC ORBIT
options.Color.Up = 'b';    %UNSTABLE PERIODIC ORBIT

%LineStyle -- SHOULD BE VALID MATLAB LINESTYLES
options.LineStyle.Ss = '-';    %STABLE STEADY STATE
options.LineStyle.Us = '-';   %UNSTABLE STEADY STATE
options.LineStyle.Sp = 'none';    %STABLE PERIODIC ORBIT
options.LineStyle.Up = 'none';    %UNSTABLE PERIODIC ORBIT

%Linewidths
options.Linewidth.Ss = 2;   %STABLE STEADY STATE
options.Linewidth.Us = 1;     %UNSTABLE STEADY STATE
options.Linewidth.Sp = 1;     %STABLE PERIODIC ORBIT
options.Linewidth.Up = 1;     %UNSTABLE PERIODIC ORBIT

%Markersize
options.Marker.Ss = 'none';    %STABLE STEADY STATE
options.Marker.Us = 'none';    %UNSTABLE STEADY STATE
options.Marker.Sp = '.';    %STABLE PERIODIC ORBIT
options.Marker.Up = 'o';    %UNSTABLE PERIODIC ORBIT

%Markersize
options.Markersize.Ss = 1;    %STABLE STEADY STATE
options.Markersize.Us = 1;    %UNSTABLE STEADY STATE
options.Markersize.Sp = 5;    %STABLE PERIODIC ORBIT
options.Markersize.Up = 5;    %UNSTABLE PERIODIC ORBIT

%MarkerFaceColor -- SHOULD BE VALID MATLAB COLORS
options.MarkerFaceColor.Ss = 'none';   %STABLE STEADY STATE
options.MarkerFaceColor.Us = 'none';   %UNSTABLE STEADY STATE
options.MarkerFaceColor.Sp = 'r';    %STABLE PERIODIC ORBIT
options.MarkerFaceColor.Up = 'none';   %UNSTABLE PERIODIC ORBIT

options.nSubSample=[1,1,1,1];

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
            options.(deblank(Names(j,:))).Ss = arg{1};
            options.(deblank(Names(j,:))).Us = arg{2};
            options.(deblank(Names(j,:))).Sp = arg{3};
            options.(deblank(Names(j,:))).Up = arg{4};
        end
        expectval = 0;
        
    end
    i = i + 1;
end

if expectval
    error(message('MATLAB:plotxppautset:NoValueForProp', arg));
end

