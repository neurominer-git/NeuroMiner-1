function k = txtmenu(xHeader,varargin)
%TXTMENU   Generate a numbered list text menu in the command window.
%   CHOICE = TXTMENU(HEADER, DEFAULT, ITEM1, ITEM2, ..., ITEMn) where all
%   arguments are strings displays
%
%   ------- HEADER -------
%  
%     [0]   DEFAULT
%      1    ITEM1
%      2    ITEM2
%      ...  ... 
%      n    ITEMn
%
%   and the user is prompted to select a menu number, which is returned as
%   CHOICE (the prompt can be changed by editing the string in the first
%   line of active code in TXTMENU). No selection (i.e. hitting return key
%   without an input) returns 0. TXTMENU will return ONLY zero or a
%   positive integer for CHOICE.
%
%   CHOICE = TXTMENU(HEADER, ITEMLIST) where ITEMLIST is a cell array of
%   the form {DEFAULT, ITEM1, ITEM2, ..., ITEMn} is also valid.
%
%   If DEFAULT is an empty matrix or string ([] or ''), the first line is
%   omitted and both 0 and no selection are illegal inputs and TXTMENU will
%   only return a positive integer.
% 
%   If the HEADER string is an empty matrix or string, the header line is
%   omitted.
%
%   ITEMs or DEFAULT may use \n (linefeed), \t (tab), etc. Linefeed is
%   useful for visually separating groups of options. Use \\ to produce a
%   backslash character and %% to produce the percent character. Call
%   SPRINTF for an item for more elaborate formatting.
%
%   Example:
%     choice = txtmenu('Main menu',...
%       {'Return / up menu\n'
%        'Change project''s name\n\t\t\t-or-\n\t\tCreate new project\n'
%        'Load project data'
%        sprintf('Save project %0.4f \n',pi)
%        'Visualize and modify data'})
%
%   Author:   Sky Sartorius
%             sky.sartorius@gmail.com
%
%   See also MENU, SPRINTF, INPUT.
prompt = 'Select a menu number: ';

if nargin < 2 || isempty(varargin{1})
    disp('No menu items to choose from.')
    k=NaN;
    return;
elseif nargin==2 && iscell(varargin{1}),
  ArgsIn = varargin{1};
else
  ArgsIn = varargin;
end

if isempty(ArgsIn{1})
    usedefault = false;
else
    usedefault = true;
end

numItems = length(ArgsIn);

if numItems == 1 && usedefault == 0;
    disp('Incomplete item list')
    k=NaN;
    return
end

while 1,
    disp(' ')
    if ~isempty(xHeader)
        disp(['------- ',sprintf(xHeader),' -------'])
    end
    disp(' ')

    if usedefault
        disp( [ '  [0]   ' sprintf(ArgsIn{1}) ] )
    end
    for n = 1 : numItems-1
        if n<10
            xgap = '   ';
        elseif n<100
            xgap = '  ';
        elseif n<1000
            xgap = ' ';
        else
            xgap = '';
        end
        disp( [ xgap '(' int2str(n) ')  ' sprintf(ArgsIn{n+1}) ] )
    end
    disp(' ')
    k = input(prompt);

    if isempty(k) || k == 0
        if usedefault
            k = 0;
        else
            k = -1;
        end
    end
        
    if  (k < 0) || (k > numItems-1) ...
        || ~strcmp(class(k),'double') ...
        || ~isreal(k) || (isnan(k)) || isinf(k) ...
        || k ~= round(k)
        disp(' ')
        disp('Selection out of range. Try again.')
        disp(' ')
    else
        disp (' ')
        disp (' ')
        return
    end
end

% Revision history:
%   V1.0    24 July 2010
%   V2.0    26 July 2010
%       Allowed for the use of \n, \t, etc. in input strings
%       Put prompt as first active line so it can be easily changed, but
%       decided not to have it be a second input for reasons of
%       compatibility with previous versions as well as MENU, and also
%       because of complicating use.
%   V3.0    4 August 2010
%       Moved default option to beginning of list. The numerical sequence
%       is maintained this way, and since the default value will often be
%       known very early in development, it saves more changing of
%       following code.
%       Allow empty header and omission of header line