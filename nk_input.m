function varargout = nk_input(varargin)
% Comprehensive graphical and command line input function
% FORMATs (given in Programmers Help)
%_______________________________________________________________________
%
% nk_input handles most forms of interactive user input for SPM.
% (File selection is handled by spm_select.m)
%
% There are five types of input: String, Evaluated, Conditions, Buttons
% and Menus:  These prompt for string input; string input which is
% evaluated to give a numerical result; selection of one item from a
% set of buttons; selection of an item from a menu.
%
% - STRING, EVALUATED & CONDITION input -
% For STRING, EVALUATED and CONDITION input types, a prompt is
% displayed adjacent to an editable text entry widget (with a lilac
% background!). Clicking in the entry widget allows editing, pressing
% <RETURN> or <ENTER> enters the result. You must enter something,
% empty answers are not accepted. A default response may be pre-specified
% in the entry widget, which will then be outlined. Clicking the border
% accepts the default value.
%
% Basic editing of the entry widget is supported *without* clicking in
% the widget, provided no other graphics widget has the focus. (If a
% widget has the focus, it is shown highlighted with a thin coloured
% line. Clicking on the window background returns the focus to the
% window, enabling keyboard accelerators.). This enables you to type
% responses to a sequence of questions without having to repeatedly
% click the mouse in the text widgets. Supported are BackSpace and
% Delete, line kill (^U). Other standard ASCII characters are appended
% to the text in the entry widget. Press <RETURN> or <ENTER> to submit
% your response.
%
% A ContextMenu is provided (in the figure background) giving access to
% relevant utilities including the facility to load input from a file
% (see spm_load.m and examples given below): Click the right button on
% the figure background.
%
% For EVALUATED input, the string submitted is evaluated in the base
% MatLab workspace (see MatLab's `eval` command) to give a numerical
% value. This permits the entry of numerics, matrices, expressions,
% functions or workspace variables. I.e.:
%   i)  - a number, vector or matrix        e.g. "[1 2 3 4]"
%                                                "[1:4]"
%                                                "1:4"
%  ii)  - an expression                     e.g. "pi^2"
%                                                "exp(-[1:36]/5.321)"
% iii)  - a function (that will be invoked) e.g. "spm_load('tmp.dat')"
%         (function must be on MATLABPATH)       "input_cov(36,5.321)"
%  iv)  - a variable from the base workspace
%                                           e.g. "tmp"
%
% The last three options provide a great deal of power: spm_load will
% load a matrix from an ASCII data file and return the results. When
% called without an argument, spm_load will pop up a file selection
% dialog. Alternatively, this facility can be gained from the
% ContextMenu. The second example assummes a custom funcion called
% input_cov has been written which expects two arguments, for example
% the following file saved as input_cov.m somewhere on the MATLABPATH
% (~/matlab, the matlab subdirectory of your home area, and the current
% directory, are on the MATLABPATH by default):
%
%       function [x] = input_cov(n,decay)
%       % data input routine - mono-exponential covariate
%       % FORMAT [x] = input_cov(n,decay)
%       % n     -  number of time points
%       % decay - decay constant
%       x = exp(-[1:n]/decay);
%
% Although this example is trivial, specifying large vectors of
% empirical data (e.g. reaction times for 72 scans) is efficient and
% reliable using this device. In the last option, a variable called tmp
% is picked up from the base workspace. To use this method, set the
% variables in the MatLab base workspace before starting an SPM
% procedure (but after starting the SPM interface). E.g.
% >> tmp=exp(-[1:36]/5.321)
%
% Occasionally a vector of a specific length will be required: This
% will be indicated in the prompt, which will start with "[#]", where
% # is the length of vector(s) required. (If a matrix is entered then
% at least one dimension should equal #.)
%
% Occasionally a specific type of number will be required. This should
% be obvious from the context. If you enter a number of the wrong type,
% you'll be alerted and asked to re-specify. The types are i) Real
% numbers; ii) Integers; iii) Whole numbers [0,1,2,3,...] & iv) Natural
% numbers [1,2,3,...]
%
% CONDITIONS type input is for getting indicator vectors. The features
% of evaluated input described above are complimented as follows:
%   v)  - a compressed list of digits 0-9   e.g. "12121212"
%  ii)  - a list of indicator characters    e.g. "abababab"
%         a-z mapped to 1-26 in alphabetical order, *except* r ("rest")
%         which is mapped to zero (case insensitive, [A:Z,a:z] only)
% ...in addition the response is checked to ensure integer condition indices.
% Occasionally a specific number of conditions will be required: This
% will be indicated in the prompt, which will end with (#), where # is
% the number of conditions required.
%
% CONTRAST type input is for getting contrast weight vectors. Enter
% contrasts as row-vectors. Contrast weight vectors will be padded with
% zeros to the correct length, and checked for validity. (Valid
% contrasts are estimable, which are those whose weights vector is in
% the row-space of the design matrix.)
%
% Errors in string evaluation for EVALUATED & CONDITION types are
% handled gracefully, the user notified, and prompted to re-enter.
%
% - BUTTON input -
% For Button input, the prompt is displayed adjacent to a small row of
% buttons. Press the approprate button. The default button (if
% available) has a dark outline. Keyboard accelerators are available
% (provided no graphics widget has the focus):  <RETURN> or <ENTER>
% selects the default button (if available). Typing the first character
% of the button label (case insensitive) "presses" that button. (If
% these Keys are not unique, then the integer keys 1,2,...  "press" the
% appropriate button.)
%
% The CommandLine variant presents a simple menu of buttons and prompts
% for a selection. Any default response is indicated, and accepted if
% an empty line is input.
%
%
% - MENU input -
% For Menu input, the prompt is displayed in a pull down menu widget.
% Using the mouse, a selection is made by pulling down the widget and
% releasing the mouse on the appropriate response. The default response
% (if set) is marked with an asterisk. Keyboard accelerators are
% available (provided no graphic widget has the focus) as follows: 'f',
% 'n' or 'd' move forward to next response down; 'b', 'p' or 'u' move
% backwards to the previous response up the list; the number keys jump
% to the appropriate response number; <RETURN> or <ENTER> slelects the
% currently displayed response. If a default is available, then
% pressing <RETURN> or <ENTER> when the prompt is displayed jumps to
% the default response.
%
% The CommandLine variant presents a simple menu and prompts for a selection.
% Any default response is indicated, and accepted if an empty line is
% input.
%
%
% - Combination BUTTON/EDIT input -
% In this usage, you will be presented with a set of buttons and an
% editable text widget. Click one of the buttons to choose that option,
% or type your response in the edit widget. Any default response will
% be shown in the edit widget. The edit widget behaves in the same way
% as with the STRING/EVALUATED input, and expects a single number.
% Keypresses edit the text widget (rather than "press" the buttons)
% (provided no other graphics widget has the focus). A default response
% can be selected with the mouse by clicking the thick border of the
% edit widget.
%
%
% - Comand line -
% If YPos is 0 or global CMDLINE is true, then the command line is used.
% Negative YPos overrides CMDLINE, ensuring the GUI is used, at
% YPos=abs(YPos). Similarly relative YPos beginning with '!'
% (E.g.YPos='!+1') ensures the GUI is used.
%
% nk_input uses the SPM 'Interactive' window, which is 'Tag'ged
% 'Interactive'. If there is no such window, then the current figure is
% used, or an 'Interactive' window created if no windows are open.
%
%-----------------------------------------------------------------------
% Programers help is contained in the main body of nk_input.m
%-----------------------------------------------------------------------
% See      : input.m     (MatLab Reference Guide)
% See also : spm_select.m   (SPM file selector dialog)
%          : nk_input.m (Input wrapper function - handles batch mode)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Andrew Holmes
% $Id: nk_input.m 816 2007-05-24 19:14:01Z karl $


%=======================================================================
% - FORMAT specifications for programers
%=======================================================================
% generic    - [p,YPos] = nk_input(Prompt,YPos,Type,...)
% string     - [p,YPos] = nk_input(Prompt,YPos,'s',DefStr)
% string+    - [p,YPos] = nk_input(Prompt,YPos,'s+',DefStr)
% evaluated  - [p,YPos] = nk_input(Prompt,YPos,'e',DefStr,n)
% - natural  - [p,YPos] = nk_input(Prompt,YPos,'n',DefStr,n,mx)
% - whole    - [p,YPos] = nk_input(Prompt,YPos,'w',DefStr,n,mx)
% - integer  - [p,YPos] = nk_input(Prompt,YPos,'i',DefStr,n)
% - real     - [p,YPos] = nk_input(Prompt,YPos,'r',DefStr,n,mm)
% condition  - [p,YPos] = nk_input(Prompt,YPos,'c',DefStr,n,m)
% contrast   - [p,YPos] = nk_input(Prompt,YPos,'x',DefStr,n,X)
% permutation- [p,YPos] = nk_input(Prompt,YPos,'p',DefStr,P,n)
% button     - [p,YPos] = nk_input(Prompt,YPos,'b',Labels,Values,DefItem)
% button/edit combo's (edit for string or typed scalar evaluated input)
%              [p,YPos] = nk_input(Prompt,YPos,'b?1',Labels,Values,DefStr,mx)
%   ...where ? in b?1 specifies edit widget type as with string & eval'd input
%            - [p,YPos] = nk_input(Prompt,YPos,'n1',DefStr,mx)
%            - [p,YPos] = nk_input(Prompt,YPos,'w1',DefStr,mx)
% button dialog 
%            - [p,YPos] = nk_input(Prompt,YPos,'bd',...
%                                Labels,Values,DefItem,Title)
% menu       - [p,YPos] = nk_input(Prompt,YPos,'m',Labels,Values,DefItem)
% display    -            nk_input(Message,YPos,'d',Label)
% display    - (GUI only) nk_input(Alert,YPos,'d!',Label)
%
% yes/no     - [p,YPos] = nk_input(Prompt,YPos,'y/n',Values,DefItem)
% buttons (shortcut) where Labels is a bar delimited string
%            - [p,YPos] = nk_input(Prompt,YPos,Labels,Values,DefItem)
%
% NB: Natural numbers are [1:Inf), Whole numbers are [0:Inf)
% 
% -- Parameters (input) --
%
% Prompt   - prompt string
%          - Defaults (missing or empty) to 'Enter an expression'
%
% YPos     - (numeric) vertical position {1 - 12}
%                  - overriden by global CMDLINE
%                  - 0 for command line
%                  - negative to force GUI
%          - (string) relative vertical position E.g. '+1'  
%                  - relative to last position used
%                  - overriden by global CMDLINE
%                  - YPos(1)=='!' forces GUI E.g. '!+1'
%                  - '_' is a shortcut for the lowest GUI position
%          - Defaults (missing or empty) to '+1'
%
% Type     - type of interrogation
%                  - 's'tring
%                  - 's+' multi-line string
%                     - p returned as cellstr (nx1 cell array of strings)
%                     - DefStr can be a cellstr or string matrix
%                  - 'e'valuated string
%                     - 'n'atural numbers
%                     - 'w'hole numbers
%                     - 'i'ntegers
%                     - 'r'eals
%                  - 'c'ondition indicator vector
%                  - 'x' - contrast entry
%                     - If n(2) or design matrix X is specified, then
%                       contrast matrices are padded with zeros to have
%                       correct length.
%                     - if design matrix X is specified, then contrasts are
%                       checked for validity (i.e. in the row-space of X)
%                       (checking handled by spm_SpUtil)
%                  - 'b'uttons
%                  - 'bd' - button dialog: Uses MatLab's questdlg
%                     - For up to three buttons
%                     - Prompt can be a cellstr with a long multiline message
%                     - CmdLine support as with 'b' type
%                  - button/edit combo's: 'be1','bn1','bw1','bi1','br1'
%                     - second letter of b?1 specifies type for edit widget
%                  - 'n1' - single natural number (buttons 1,2,... & edit)
%                  - 'w1' - single whole number   (buttons 0,1,... & edit)
%                  - 'm'enu pulldown
%                  - 'y/n' : Yes or No buttons
%                                    (See shortcuts below)
%                  - bar delimited string : buttons with these labels
%                                    (See shortcuts below)
%          - Defaults (missing or empty) to 'e'
%
% DefStr   - Default string to be placed in entry widget for string and
%            evaluated types
%          - Defaults to ''
%
% n ('e', 'c' & 'p' types)
%          - Size of matrix requred
%          - NaN for 'e' type implies no checking - returns input as evaluated
%          - length of n(:) specifies dimension - elements specify size
%          - Inf implies no restriction
%          - Scalar n expanded to [n,1] (i.e. a column vector)
%            (except 'x' contrast type when it's [n,np] for np
%          - E.g: [n,1] & [1,n] (scalar n) prompt for an n-vector,
%                         returned as column or row vector respectively
%                 [1,Inf] & [Inf,1] prompt for a single vector,
%                         returned as column or row vector respectively
%                 [n,Inf] & [Inf,n] prompts for any number of n-vectors,
%                         returned with row/column dimension n respectively.
%                 [a,b] prompts for an 2D matrix with row dimension a and
%                         column dimension b
%                 [a,Inf,b] prompt for a 3D matrix with row dimension a,
%                         page dimension b, and any column dimension.
%          - 'c' type can only deal with single vectors
%          - NaN for 'c' type treated as Inf
%          - Defaults (missing or empty) to NaN
%
% n ('x'type)
%          - Number of contrasts required by 'x' type (n(1))
%            ( n(2) can be used specify length of contrast vectors if )
%            ( a design matrix isn't passed                           )
%          - Defaults (missing or empty) to 1 - vector contrast
%
% mx ('n', 'w', 'n1', 'w1', 'bn1' & 'bw1' types)
%          - Maximum value (inclusive)
%
% mm ('r' type)
%          - Maximum and minimum values (inclusive)
%
% m        - Number of unique conditions required by 'c' type
%          - Inf implies no restriction
%          - Defaults (missing or empty) to Inf - no restriction
%
% P        - set (vector) of numbers of which a permutation is required
%
% X        - Design matrix for contrast checking in 'x' type
%          - Can be either a straight matrix or a space structure (see spm_sp)
%          - Column dimension of design matrix specifies length of contrast
%            vectors (overriding n(2) is specified).
%
% Title    - Title for questdlg in 'bd' type
%
% Labels   - Labels for button and menu types.
%                  - string matrix, one label per row
%                  - bar delimited string
%                            E.g. 'AnCova|Scaling|None'
%
% Values   - Return values corresponding to Labels for button and menu types
%          - j-th row is returned if button / menu item j is selected
%            (row vectors are transposed)
%          - Defaults (missing or empty) to - (button) Labels
%                                           - ( menu ) menu item numbers
%
% DefItem  - Default item number, for button and menu types.
%
% -- Parameters (output) --
% p        - results
% YPos     - Optional second output argument returns GUI position just used
%
%-----------------------------------------------------------------------
% WINDOWS:
%
% nk_input uses the SPM 'Interactive' 'Tag'ged window. If this isn't
% available and no figures are open, an 'Interactive' SPM window is
% created (`spm('CreateIntWin')`). If figures are available, then the
% current figure is used *unless* it is 'Tag'ged.
%
%-----------------------------------------------------------------------
% SHORTCUTS:
%
% Buttons SHORTCUT - If the Type parameter is a bar delimited string, then
% the Type is taken as 'b' with the specified labels, and the next parameter
% (if specified) is taken for the Values.
%
% Yes/No question shortcut - p = nk_input(Prompt,YPos,'y/n') expands
% to p = nk_input(Prompt,YPos,'b','yes|no',...), enabling easy use of
% nk_input for yes/no dialogue. Values defaults to 'yn', so 'y' or 'n'
% is returned as appropriate.
%
%-----------------------------------------------------------------------
% EXAMPLES:
%            ( Specified YPos is overriden if global CMDLINE is )
%            ( true, when the command line versions are used.   )
%
%       p = nk_input
%               Command line input of an evaluated string, default prompt.
%       p = nk_input('Enter a value',1)
%               Evaluated string input, prompted by 'Enter a value', in
%               position 1 of the dialog figure.
%       p = nk_input(str,'+1','e',0.001)
%               Evaluated string input, prompted by contents of string str,
%               in next position of the dialog figure.
%               Default value of 0.001 offered.
%       p = nk_input(str,2,'e',[],5)
%               Evaluated string input, prompted by contents of string str,
%               in second position of the dialog figure.
%               Vector of length 5 required - returned as column vector
%       p = nk_input(str,2,'e',[],[Inf,5])
%               ...as above, but can enter multiple 5-vectors in a matrix,
%               returned with 5-vectors in rows
%       p = nk_input(str,0,'c','ababab')
%               Condition string input, prompted by contents of string str
%               Uses command line interface.
%               Default string of 'ababab' offered.
%       p = nk_input(str,0,'c','010101')
%               As above, but default string of '010101' offered.
%       [p,YPos] = nk_input(str,'0','s','Image')
%               String input, same position as last used, prompted by str,
%               default of 'Image' offered. YPos returns GUI position used.
%       p = nk_input(str,'-1','y/n')
%               Yes/No buttons for question with prompt str, in position one 
%               before the last used Returns 'y' or 'n'.
%       p = nk_input(str,'-1','y/n',[1,0],2)
%               As above, but returns 1 for yes response, 0 for no,
%               with 'no' as the default response
%       p = nk_input(str,4,'AnCova|Scaling')
%               Presents two buttons labelled 'AnCova' & 'Scaling', with 
%               prompt str, in position 4 of the dialog figure. Returns the 
%               string on the depresed button, where buttons can be pressed 
%               with the mouse or by the respective keyboard accelerators
%               'a' & 's' (or 'A' & 'S').
%       p = nk_input(str,-4,'b','AnCova|Scaling',[],2)
%               As above, but makes "Scaling" the default response, and
%               overrides global CMDLINE
%       p = nk_input(str,0,'b','AnCova|Scaling|None',[1,2,3])
%               Prompts for [A]ncova / [S]caling / [N]one in MatLab command
%               window, returns 1, 2, or 3 according to the first character
%               of the entered string as one of 'a', 's', or 'n' (case 
%               insensitive).
%       p = nk_input(str,1,'b','AnCova',1)
%		Since there's only one button, this just displays the response
%		in GUI position 1 (or on the command line if global CMDLINE
%		is true), and returns 1.
%	p = nk_input(str,'+0','br1','None|Mask',[-Inf,NaN],0.8)
%               Presents two buttons labelled "None" & "Mask" (which return
%               -Inf & NaN if clicked), together with an editable text widget
%               for entry of a single real number. The default of 0.8 is
%               initially presented in the edit window, and can be selected by
%               pressing return.
%		Uses the previous GUI position, unless global CMDLINE is true,
%               in which case a command-line equivalent is used.
%	p = nk_input(str,'+0','w1')
%		Prompts for a single whole number using a combination of
%		buttons and edit widget, using the previous GUI position,
%		or the command line if global CMDLINE is true.
%       p = nk_input(str,'!0','m','Single Subject|Multi Subject|Multi Study')
%               Prints the prompt str in a pull down menu containing items
%               'Single Subject', 'Multi Subject' & 'Multi Study'. When OK is
%               clicked p is returned as the index of the  choice, 1,2, or 3 
%               respectively. Uses last used position in GUI, irrespective of 
%               global CMDLINE
%       p = nk_input(str,5,'m',...
%               'Single Subject|Multi Subject|Multi Study',...
%               ['SS';'MS';'SP'],2)
%               As above, but returns strings 'SS', 'MS', or 'SP' according to
%               the respective choice, with 'MS; as the default response.
%       p = nk_input(str,0,'m',...
%               'Single Subject|Multi Subject|Multi Study',...
%               ['SS';'MS';'SP'],2)
%               As above, but the menu is presented in the command window
%               as a numbered list.
%       nk_input('AnCova, GrandMean scaling',0,'d')
%               Displays message in a box in the MatLab command window
%       [null,YPos]=nk_input('Session 1','+1','d!','fMRI')
%		Displays 'fMRI: Session 1' in next GUI position of the
%               'Interactive' window. If CMDLINE is 1, then nothing is done.
%               Position used is returned in YPos.
%
%-----------------------------------------------------------------------
% FORMAT h = nk_input(Prompt,YPos,'m!',Labels,cb,UD,XCB);
% GUI PullDown menu utility - creates a pulldown menu in the Interactive window
% FORMAT H = nk_input(Prompt,YPos,'b!',Labels,cb,UD,XCB);
% GUI Buttons utility - creates GUI buttons in the Interactive window
%
% Prompt, YPos, Labels - as with 'm'enu/'b'utton types
% cb  - CallBack string
% UD  - UserData
% XCB - Extended CallBack handling - allows different CallBack for each item,
%       and use of UD in CallBack strings. [Defaults to 1 for PullDown type
%       when multiple CallBacks specified, 0 o/w.]
% H   - Handle of 'PullDown' uicontrol / 'Button's
%
% In "normal" mode (when XCB is false), this is essentially a utility
% to create a PullDown menu widget or set of buttons in the SPM
% 'Interactive' figure, using positioning and Label definition
% conveniences of the nk_input 'm'enu & 'b'utton types. If Prompt is
% not empty, then the PullDown/Buttons appears on the right, with the
% Prompt on the left, otherwise the PullDown/Buttons use the whole
% width of the Interactive figure. The PopUp's CallBack string is
% specified in cb, and [optional] UserData may be passed as UD.
%
% For buttons, a separate callback can be specified for each button, by
% passing the callbacks corresponding to the Labels as rows of a
% cellstr or string matrix.
%
% This "different CallBacks" facility can also be extended to the
% PullDown type, using the "extended callback" mode (when XCB is
% true).  % In addition, in "extended callback", you can use UD to
% refer to the UserData argument in the CallBack strings. (What happens
% is this: The cb & UD are stored as fields in the PopUp's UserData
% structure, and the PopUp's callback is set to nk_input('!m_cb'),
% which reads UD into the functions workspace and eval's the
% appropriate CallBack string.  Note that this means that base
% workspace variables are inaccessible (put what you need in UD), and
% that any return arguments from CallBack functions are not passed back
% to the base workspace).
% 
%
%-----------------------------------------------------------------------
% UTILITY FUNCTIONS:
%
% FORMAT colour = nk_input('!Colour')
% Returns colour for input widgets, as specified in COLOUR parameter at 
% start of code.
% colour  - [r,g,b] colour triple
%
% FORMAT [iCond,msg] = nk_input('!iCond',str,n,m)
% Parser for special 'c'ondition type: Handles digit strings and
% strings of indicator chars.
% str     - input string
% n       - length of condition vector required	[defaut Inf - no restriction]
% m       - number of conditions required	[default Inf - no restrictions]
% iCond   - Integer condition indicator vector
% msg     - status message
%
% FORMAT hM = nk_input('!InptConMen',Finter,H)
% Sets a basic Input ContextMenu for the figure
% Finter - figure to set menu in
% H      - handles of objects to delete on "crash out" option
% hM     - handle of UIContextMenu
%
% FORMAT [CmdLine,YPos] = nk_input('!CmdLine',YPos)
% Sorts out whether to use CmdLine or not & canonicalises YPos
% CmdLine - Binary flag
% YPos    - Position index
%
% FORMAT Finter = nk_input('!GetWin',F)
% Locates (or creates) figure to work in
% F       - Interactive Figure, defaults to 'Interactive'
% Finter  - Handle of figure to use
%
% FORMAT [PLoc,cF] = nk_input('!PointerJump',RRec,F,XDisp)
% Raise window & jump pointer over question
% RRec  - Response rectangle of current question
% F     - Interactive Figure, Defaults to 'Interactive'
% XDisp - X-displacement of cursor relative to RRec
% PLoc  - Pointer location before jumping
% cF    - Current figure before making F current.
%
% FORMAT [PLoc,cF] = nk_input('!PointerJumpBack',PLoc,cF)
% Replace pointer and reset CurrentFigure back
% PLoc  - Pointer location before jumping
% cF    - Previous current figure
%
% FORMAT nk_input('!PrntPrmpt',Prompt,TipStr,Title)
% Print prompt for CmdLine questioning
% Prompt - prompt string, callstr, or string matrix
% TipStr - tip string
% Title  - title string
%
% FORMAT [Frec,QRec,PRec,RRec] = nk_input('!InputRects',YPos,rec,F)
% Returns rectangles (pixels) used in GUI
% YPos  - Position index
% rec   - Rectangle specifier: String, one of 'Frec','QRec','PRec','RRec'
%         Defaults to '', which returns them all.
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
% FRec  - Position of interactive window
% QRec  - Position of entire question
% PRec  - Position of prompt
% RRec  - Position of response
%
% FORMAT nk_input('!DeleteInputObj',F)
% Deltes input objects (only) from figure F
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
%
% FORMAT [CPos,hCPos] = nk_input('!CurrentPos',F)
% Returns currently used GUI question positions & their handles
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
% CPos  - Vector of position indices
% hCPos - (n x CPos) matrix of object handles
%
% FORMAT h = nk_input('!FindInputObj',F)
% Returns handles of input GUI objects in figure F
% F - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
% h - vector of object handles
%
% FORMAT [NPos,CPos,hCPos] = nk_input('!NextPos',YPos,F,CmdLine)
% Returns next position index, specified by YPos
% YPos    - Absolute (integer) or relative (string) position index
%           Defaults to '+1'
% F       - Interactive Figure, defaults to spm_figure('FindWin','Interactive')
% CmdLine - Command line? Defaults to nk_input('!CmdLine',YPos)
% NPos    - Next position index
% CPos & hCPos - as for !CurrentPos
%
% FORMAT NPos = nk_input('!SetNextPos',YPos,F,CmdLine)
% Sets up for input at next position index, specified by YPos. This utility
% function can be used stand-alone to implicitly set the next position
% by clearing positions NPos and greater.
% YPos    - Absolute (integer) or relative (string) position index
%           Defaults to '+1'
% F       - Interactive Figure, defaults to spm_figure('FindWin','Interactive')
% CmdLine - Command line? Defaults to nk_input('!CmdLine',YPos)
% NPos    - Next position index
%
% FORMAT MPos = nk_input('!MaxPos',F,FRec3)
% Returns maximum position index for figure F
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
%	  Not required if FRec3 is specified
% FRec3 - Length of interactive figure in pixels
% 
% FORMAT nk_input('!EditableKeyPressFcn',h,ch)
% KeyPress callback for GUI string / eval input
%
% FORMAT nk_input('!ButtonKeyPressFcn',h,Keys,DefItem,ch)
% KeyPress callback for GUI buttons
%
% FORMAT nk_input('!PullDownKeyPressFcn',h,ch,DefItem)
% KeyPress callback for GUI pulldown menus
%
% FORMAT nk_input('!m_cb')
% Extended CallBack handler for 'p' PullDown utility type
%
% FORMAT nk_input('!dScroll',h,str)
% Scroll text string in object h
% h      - handle of text object
% Prompt - Text to scroll (Defaults to 'UserData' of h)
%
%-----------------------------------------------------------------------
% SUBFUNCTIONS:
%
% FORMAT [Keys,Labs] = sf_labkeys(Labels)
% Make unique character keys for the Labels, ignoring case.
% Used with 'b'utton types.
%
% FORMAT [p,msg] = sf_eEval(str,Type,n,m)
% Common code for evaluating various input types.
%
% FORMAT str = sf_SzStr(n,l)
% Common code to construct prompt strings for pre-specified vector/matrix sizes
%
% FORMAT [p,msg] = sf_SzChk(p,n,msg)
% Common code to check (& canonicalise) sizes of input vectors/matrices
%
%_______________________________________________________________________
% @(#)nk_input.m	2.8 Andrew Holmes 03/03/04
global NMinfo

%-Parameters
%=======================================================================
%COLOUR   = get(0,'defaultUicontrolBackgroundColor');
PJump    = 1;		%-Jumping of pointer to question?
TTips    = 1;		%-Use ToolTipStrings? (which can be annoying!)
ConCrash = 1;		%-Add "crash out" option to 'Interactive'fig.ContextMenu
CmdLine		=true;

%-Condition arguments
%=======================================================================
if nargin<1|isempty(varargin{1}), Prompt=''; else, Prompt=varargin{1}; end

if ~isempty(Prompt) & ischar(Prompt) & Prompt(1)=='!'
	%-Utility functions have Prompt string starting with '!'
	Type = Prompt;
else			%-Should be an input request: get Type & YPos
	if nargin<3|isempty(varargin{3}), Type='e';  else, Type=varargin{3}; end
	if any(Type=='|'), Type='b|'; end
	if nargin<2|isempty(varargin{2}), YPos='+1'; else, YPos=varargin{2}; end

	%[CmdLine,YPos] = nk_input('!CmdLine',YPos);
	
%  	if ~CmdLine	%-Setup for GUI use
%  		%-Locate (or create) figure to work in
%  		Finter = nk_input('!GetWin');
%  	
%  		%-Find out which Y-position to use, setup for use
%  		YPos = nk_input('!SetNextPos',YPos,Finter,CmdLine);
%  		
%  		%-Determine position of objects
%  		[FRec,QRec,PRec,RRec]=nk_input('!InputRects',YPos,'',Finter);
%  	end
end


switch lower(Type)
case {'s','s+','sq','e','n','w','i','r','c','x','p'}  %-String and evaluated input
%=======================================================================
%-Condition arguments
if nargin<7|isempty(varargin{7}), l=[]; else, l=varargin{7}; end
if nargin<6|isempty(varargin{6}), m=[]; else, m=varargin{6}; end
if nargin<5|isempty(varargin{5}), n=[]; else, n=varargin{5}; end
if nargin<4, DefStr=''; else, DefStr=varargin{4}; end
if strcmp(lower(Type),'s+')
	%-DefStr should be a cellstr for 's+' type.
	if isempty(DefStr), DefStr = {};
		else, DefStr = cellstr(DefStr); end
	DefStr = DefStr(:);
else
	%-DefStr needs to be a string
	if ~ischar(DefStr), DefStr=num2str(DefStr); end
	DefStr = DefStr(:)';
end

strM='';
switch lower(Type)			%-Type specific defaults/setup
case 's', TTstr='enter string';
case 's+',TTstr='enter string - multi-line';
case 'sq',TTstr='enter string or enter to return';
case 'e', TTstr='enter expression to evaluate';
case 'n', TTstr='enter expression - natural number(s)';
	if ~isempty(m), strM=sprintf(' (in [1,%d])',m); TTstr=[TTstr,strM]; end
case 'w', TTstr='enter expression - whole number(s)';
	if ~isempty(m) && ~isempty(l), 
        strM=sprintf(' (in [%d,%d])',l,m); TTstr=[TTstr,strM]; 
    elseif ~isempty(m)
        strM=sprintf(' (in [0,%d])',m); TTstr=[TTstr,strM]; 
    end
case 'i', TTstr='enter expression - integer(s)';
case 'r', TTstr='enter expression - real number(s)';
	if ~isempty(m), TTstr=[TTstr,sprintf(' in [%g,%g]',min(m),max(m))]; end
case 'c', TTstr='enter indicator vector e.g. 0101...  or abab...';
	if ~isempty(m) & isfinite(m), strM=sprintf(' (%d)',m); end
case 'x', TTstr='enter contrast matrix';
case 'p',
	if isempty(n), error('permutation of what?'), else, P=n(:)'; end
	if isempty(m), n = [1,length(P)]; end
	m = P;
	if ~length(setxor(m,[1:max(m)]))
		TTstr=['enter permutation of [1:',num2str(max(m)),']'];
	else
		TTstr=['enter permutation of [',num2str(m),']'];
	end
otherwise, TTstr='enter expression'; end

strN = sf_SzStr(n);


%  if CmdLine                                   %-Use CmdLine to get answer
%-----------------------------------------------------------------------
	nk_input('!PrntPrmpt',[Prompt,strN,strM],TTstr)

	%-Do Eval Types in Base workspace, catch errors
	switch lower(Type), case 's'
		if ~isempty(DefStr)
			Prompt=[Prompt,' (Default: ',DefStr,' )'];
		end
		str = input([Prompt,' : '],'s');
		if isempty(str), str=DefStr; end

		while isempty(str)
			spm('Beep')
			fprintf('! %s : enter something!\n',mfilename)
			str = input([Prompt,' : '],'s');
			if isempty(str), str=DefStr; end
		end
		p = str; msg = '';

	case 's+'
		fprintf(['Multi-line input: Type ''.'' on a line',...
			' of its own to terminate input.\n'])
		if ~isempty(DefStr)
			fprintf('Default : (press return to accept)\n')
			fprintf('        : %s\n',DefStr{:})
		end
		fprintf('\n')

		str = input('l001 : ','s');
		while (isempty(str) | strcmp(str,'.')) & isempty(DefStr)
			spm('Beep')
			fprintf('! %s : enter something!\n',mfilename)
			str = input('l001 : ','s');
		end

		if isempty(str)
			%-Accept default
			p = DefStr;
		else
			%-Got some input, allow entry of additional lines
			p = {str};
			str = input(sprintf('l%03u : ',length(p)+1),'s');
			while ~strcmp(str,'.')
				p = [p;{str}];
				str = input(sprintf('l%03u : ',length(p)+1),'s');
			end
		end
		msg = '';
        
        case 'sq'
            if ~isempty(DefStr)
			Prompt=[Prompt,' (Default: ',DefStr,' )'];
            end
            str = input([Prompt,' : '],'s');
            if isempty(str), str=DefStr; end
            p = str; msg = '';

	otherwise
		if ~isempty(DefStr)
			Prompt=[Prompt,' (Default: ',DefStr,' )'];
		end
		str = input([Prompt,' : '],'s');
		if isempty(str), str=DefStr; end
		[p,msg] = sf_eEval(str,Type,n,m,l);

		while ischar(p)
			spm('Beep'), fprintf('! %s : %s\n',mfilename,msg)
			str = input([Prompt,' : '],'s');
			if isempty(str), str=DefStr; end
			[p,msg] = sf_eEval(str,Type,n,m,l);
		end
	end
	if ~isempty(msg), fprintf('\t%s\n',msg), end

%  else                                             %-Use GUI to get answer
%  %-----------------------------------------------------------------------
%  
%  	%-Create text and edit control objects
%  	%---------------------------------------------------------------
%  	hPrmpt = uicontrol(Finter,'Style','Text',...
%  		'String',[strN,Prompt,strM],...
%  		'Tag',['GUIinput_',int2str(YPos)],...
%  		'UserData','',...
%          'BackgroundColor',COLOUR,...
%  		'HorizontalAlignment','Right',...
%  		'Position',PRec);
%  
%  
%  
%  	%-Default button surrounding edit widget (if a DefStr given)
%  	%-Callback sets hPrmpt UserData, and EditWidget string, to DefStr
%  	% (Buttons UserData holds handles [hPrmpt,hEditWidget], set later)
%  	cb = ['set(get(gcbo,''UserData'')*[1;0],''UserData'',',...
%  			'get(gcbo,''String'')),',...
%  		'set(get(gcbo,''UserData'')*[0;1],''String'',',...
%  			'get(gcbo,''String''))'];
%  	if ~isempty(DefStr)
%  		if iscellstr(DefStr), str=[DefStr{1},'...'];
%  			else, str=DefStr; end
%  		hDef = uicontrol(Finter,'Style','PushButton',...
%  			'String',DefStr,...
%  			'ToolTipString',...
%  				['Click on border to accept default: ' str],...
%  			'Tag',['GUIinput_',int2str(YPos)],...
%  			'UserData',[],...
%              'BackgroundColor',COLOUR,...
%  			'CallBack',cb,...
%  			'Position',RRec+[-2,-2,+4,+4]);
%  	else
%  		hDef = [];
%  	end
%  
%  	%-Edit widget: Callback puts string into hPrompts UserData
%  	cb = 'set(get(gcbo,''UserData''),''UserData'',get(gcbo,''String''))';
%  	h = uicontrol(Finter,'Style','Edit',...
%  		'String',DefStr,...
%  		'Max',strcmp(lower(Type),'s+')+1,...
%  		'Tag',['GUIinput_',int2str(YPos)],...
%  		'UserData',hPrmpt,...
%  		'CallBack',cb,...
%  		'Horizontalalignment','Left',...
%  		'BackgroundColor','w',...
%  		'Position',RRec);
%  	set(hDef,'UserData',[hPrmpt,h])
%  	uifocus(h);
%  	if TTips, set(h,'ToolTipString',TTstr), end
%  
%  	%-Figure ContextMenu for shortcuts
%  	hM = nk_input('!InptConMen',Finter,[hPrmpt,hDef,h]);
%  	cb = [	'set(get(gcbo,''UserData''),''String'',',...
%  			'[''spm_load('''''',spm_select(1),'''''')'']), ',...
%  		'set(get(get(gcbo,''UserData''),''UserData''),''UserData'',',...
%  			'get(get(gcbo,''UserData''),''String''))'];
%  	uimenu(hM,'Label','load from text file','Separator','on',...
%  		'CallBack',cb,'UserData',h)
%  
%  	%-Bring window to fore & jump pointer to edit widget
%  	[PLoc,cF] = nk_input('!PointerJump',RRec,Finter);
%  
%  	%-Setup FigureKeyPressFcn for editing of entry widget without clicking
%  	set(Finter,'KeyPressFcn',[...
%  	    'nk_input(''!EditableKeyPressFcn'',',...
%  	    'findobj(gcf,''Tag'',''GUIinput_',int2str(YPos),''',',...
%  	    	'''Style'',''edit''),',...
%  	    'get(gcbf,''CurrentCharacter''))'])
%  
%  
%  	%-Wait for edit, do eval Types in Base workspace, catch errors
%  	%---------------------------------------------------------------
%  	waitfor(hPrmpt,'UserData')
%  	if ~ishandle(hPrmpt), error(['Input window cleared whilst waiting ',...
%  		'for response: Bailing out!']), end
%  	str = get(hPrmpt,'UserData');
%  	switch lower(Type), case 's'
%  		p = str; msg = '';
%  	case 's+'
%  		p = cellstr(str); msg = '';
%  	otherwise
%  		[p,msg] = sf_eEval(str,Type,n,m);
%  		while ischar(p)
%  			set(h,'Style','Text',...
%  				'String',msg,'HorizontalAlignment','Center',...
%  				'ForegroundColor','r')
%  			spm('Beep'), pause(2)
%  			set(h,'Style','Edit',...
%  				'String',str,...
%  				'HorizontalAlignment','Left',...
%  				'ForegroundColor','k')
%  			%set(hPrmpt,'UserData','');
%  			waitfor(hPrmpt,'UserData')
%  			if ~ishandle(hPrmpt), error(['Input window cleared ',...
%  				'whilst waiting for response: Bailing out!']),end
%  			str = get(hPrmpt,'UserData');
%  			[p,msg] = sf_eEval(str,Type,n,m);
%  		end
%  	end
%  
%  	%-Fix edit window, clean up, reposition pointer, set CurrentFig back
%  	delete([hM,hDef]), set(Finter,'KeyPressFcn','')
%  	set(h,'Style','Text','HorizontalAlignment','Center',...
%  		'ToolTipString',msg,...
%  		'BackgroundColor',COLOUR)
%  	nk_input('!PointerJumpBack',PLoc,cF)
%  	drawnow
%  
%  end % (if CmdLine)

%-Return response
%-----------------------------------------------------------------------
varargout = {p,YPos};


case {'b','bd','b|','y/n','be1','bn1','bw1','bi1','br1',...
	'-n1','n1','-w1','w1','m','mq'}             %-'b'utton & 'm'enu Types
%=======================================================================
%-Condition arguments
switch lower(Type), case {'b','be1','bi1','br1','m', 'mq'}
	m = []; Title = '';
	if nargin<6, DefItem=[];  else, DefItem=varargin{6}; end
	if nargin<5, Values=[];   else, Values =varargin{5}; end
	if nargin<4, Labels='';   else, Labels =varargin{4}; end
case 'bd'
	if nargin<7, Title='';    else, Title  =varargin{7}; end
	if nargin<6, DefItem=[];  else, DefItem=varargin{6}; end
	if nargin<5, Values=[];   else, Values =varargin{5}; end
	if nargin<4, Labels='';   else, Labels =varargin{4}; end
case 'y/n'
	Title = '';
	if nargin<5, DefItem=[];  else, DefItem=varargin{5}; end
	if nargin<4, Values=[];   else, Values =varargin{4}; end
	if isempty(Values), Values='yn'; end
	Labels = {'yes','no'};
case 'b|'
	Title = '';
	if nargin<5, DefItem=[];  else, DefItem=varargin{5}; end
	if nargin<4, Values=[];   else, Values =varargin{4}; end
	Labels = varargin{3};
case 'bn1'
	if nargin<7, m=[];        else, m=varargin{7};       end
	if nargin<6, DefItem=[];  else, DefItem=varargin{6}; end
	if nargin<5, Values=[];   else, Values =varargin{5}; end
	if nargin<4, Labels=[1:5]'; Values=[1:5]; Type='-n1';
		else, Labels=varargin{4}; end
case 'bw1'
	if nargin<7, m=[];        else, m=varargin{7};       end
	if nargin<6, DefItem=[];  else, DefItem=varargin{6}; end
	if nargin<5, Values=[];   else, Values =varargin{5}; end
	if nargin<4, Labels=[0:4]'; Values=[0:4]; Type='-w1';
		else, Labels=varargin{4}; end
case {'-n1','n1','-w1','w1'}
	if nargin<5, m=[];        else, m=varargin{5};       end
	if nargin<4, DefItem=[];  else, DefItem=varargin{4}; end
	switch lower(Type)
	case {'n1','-n1'}, Labels=[1:min([5,m])]'; Values=Labels'; Type='-n1';
	case {'w1','-w1'}, Labels=[0:min([4,m])]'; Values=Labels'; Type='-w1';
	end
end

%-Check some labels were specified
if isempty(Labels), error('No Labels specified'), end

if iscellstr(Labels), Labels=char(Labels); end

%-Convert Labels "option" string to string matrix if required
if ischar(Labels) & any(Labels=='|')
	OptStr=Labels;
	BarPos=find([OptStr=='|',1]);
	Labels=OptStr(1:BarPos(1)-1);
	for Bar = 2:sum(OptStr=='|')+1
		Labels=strvcat(Labels,OptStr(BarPos(Bar-1)+1:BarPos(Bar)-1));
	end
end

%-Set default Values for the Labels
if isempty(Values)
	if strcmp(lower(Type),'m')
		Values=[1:size(Labels,1)]';
	else
		Values=Labels;
	end
else
	%-Make sure Values are in rows
	if size(Labels,1)>1 & size(Values,1)==1, Values = Values'; end
	%-Check numbers of Labels and Values match
	if (size(Labels,1)~=size(Values,1))
		error('Labels & Values incompatible sizes'), end
end

%-Numeric Labels to strings
if isnumeric(Labels)
	tmp = Labels; Labels = cell(size(tmp,1),1);
	for i=1:prod(size(tmp)), Labels{i}=num2str(tmp(i,:)); end
	Labels=char(Labels);
end

switch lower(Type), 
    
    case {'b','bd','b|','y/n'}    %-Process button types
    %=======================================================================
	
	%-Make unique character keys for the Labels, sort DefItem
	%---------------------------------------------------------------
	nLabels     = size(Labels,1);
	[Keys,Labs] = sf_labkeys(Labels);

	if ~isempty(DefItem) & any(DefItem==[1:nLabels])
		DefKey = Keys(DefItem);
	else
		DefItem = 0;
		DefKey  = '';
	end

    %-Display question prompt
    nk_input('!PrntPrmpt',Prompt,'',Title)
    %-Build prompt
    %-------------------------------------------------------
    if ~isempty(Labs) 
        Prmpt = ['[',Keys(1),']',deblank(Labs(1,:)),' '];
        for i = 2:nLabels
            Prmpt=[Prmpt,'/ [',Keys(i),']',deblank(Labs(i,:)),' '];
        end
    else
        Prmpt = ['[',Keys(1),'] '];
        for i = 2:nLabels, Prmpt=[Prmpt,'/ [',Keys(i),'] ']; end
    end
    if DefItem
        Prmpt = [Prmpt,...
            ' (Default: ',deblank(Labels(DefItem,:)),')'];
    end

    %-Ask for user response
    %-------------------------------------------------------
    if nLabels==1
        %-Only one choice - auto-pick & display
        k = 1; fprintf('%s: %s\t(only option)',Prmpt,Labels)
    else
        str = input([Prmpt,'? '],'s');
        if isempty(str), str=DefKey; end
        while isempty(str) | ~any(lower(Keys)==lower(str(1)))
            if ~isempty(str),fprintf('%c\t!Out of range\n',7),end
            str = input([Prmpt,'? '],'s');
            if isempty(str), str=DefKey; end
        end
        k = find(lower(Keys)==lower(str(1)));
    end
    fprintf('\n')

    p = Values(k,:); if ischar(p), p=deblank(p); end

case {'be1','bn1','bw1','bi1','br1','-n1','-w1'}
                                      %-Process button/entry combo types
    %=======================================================================
    if ischar(DefItem), DefStr=DefItem; else, DefStr=num2str(DefItem); end
    if isempty(m), strM=''; else, strM=sprintf(' (<=%d)',m); end


	%-Process default item
	%---------------------------------------------------------------
	if ~isempty(DefItem)
		[DefVal,msg] = sf_eEval(DefStr,Type(2),1);
		if ischar(DefVal), error(['Invalid DefItem: ',msg]), end
		Labels  = strvcat(Labels,DefStr);
		Values  = [Values;DefVal];
		DefItem = size(Labels,1);
	end

	%-Add option to specify...
	Labels = strvcat(Labels,'specify...');

	%-Process options
	nLabels     = size(Labels,1);
	[Keys,Labs] = sf_labkeys(Labels);

	if ~isempty(DefItem), DefKey = Keys(DefItem); else, DefKey = ''; end

	%-Print banner prompt
	%---------------------------------------------------------------
	nk_input('!PrntPrmpt',Prompt)		%-Display question prompt


	if Type(1)=='-'		%-No special buttons - go straight to input

		k = size(Labels,1);

	else			%-Offer buttons, default or "specify..."
	
		%-Build prompt
		%-------------------------------------------------------
		if ~isempty(Labs) 
		    Prmpt = ['[',Keys(1),']',deblank(Labs(1,:)),' '];
		    for i = 2:nLabels
		        Prmpt=[Prmpt,'/ [',Keys(i),']',deblank(Labs(i,:)),' '];
		    end
		else
		    Prmpt = ['[',Keys(1),'] '];
		    for i = 2:nLabels, Prmpt=[Prmpt,'/ [',Keys(i),'] ']; end
		end
		if DefItem, Prmpt = [Prmpt,...
			' (Default: ',deblank(Labels(DefItem,:)),')']; end
	
		%-Ask for user response
		%-------------------------------------------------------
		if nLabels==1
			%-Only one choice - auto-pick & display
			k = 1; fprintf('%s: %s\t(only option)',Prmpt,Labels)
		else
			str = input([Prmpt,'? '],'s');
			if isempty(str), str=DefKey; end
			while isempty(str) | ~any(lower(Keys)==lower(str(1)))
			    if ~isempty(str),fprintf('%c\t!Invalid response\n',7),end
			    str = input([Prmpt,'? '],'s');
			    if isempty(str), str=DefKey; end
			end
			k = find(lower(Keys)==lower(str(1)));
		end
		fprintf('\n')

	end

	%-Process response: prompt for value if "specify..." option chosen
	%===============================================================
	if k<size(Labels,1)
		p = Values(k,:); if ischar(p), p=deblank(p); end
	else

		%-"specify option chosen: ask user to specify
		%-------------------------------------------------------
		switch lower(Type(2))
                case 's', tstr=' string';        case 'e', tstr='n expression';
                case 'n', tstr=' natural number';case 'w', tstr=' whole number';
                case 'i', tstr='n integer';      case 'r', tstr=' real number';
                otherwise, tstr=''; 
        end
		
		Prompt = sprintf('%s (a%s%s)',Prompt,tstr,strM);
		if ~isempty(DefStr)
			Prompt=sprintf('%s\b, default %s)',Prompt,DefStr); end
		str = input([Prompt,' : '],'s');
		if isempty(str), str=DefStr; end
	
		%-Eval in Base workspace, catch errors
		[p,msg] = sf_eEval(str,Type(2),1,m);
     
		while ischar(p)
			spm('Beep'), fprintf('! %s : %s\n',mfilename,msg)
			str = input([Prompt,' : '],'s');
			if isempty(str), str=DefStr; end
			[p,msg] = sf_eEval(str,Type(2),1,m);
		end
	end

case {'m','mq'}                                             %-Process menu type
%=======================================================================
	
    Labels = nk_FormatNicely(Labels);
    LabelWidth = size(Labels,2); 
    PromptWidth = size(Prompt,2);
    if PromptWidth>LabelWidth, LabelWidth=PromptWidth; end
    nk_input('!PrntPrmpt',Prompt,[],[],LabelWidth)
    nLabels = size(Labels,1);
    switch lower(Type)
        case 'm'
            for i = 1:nLabels, fprintf('\t%2d | %s\n',i,Labels(i,:)), end
            Prmpt = ['Menu choice (1-',int2str(nLabels),')'];
            if DefItem
                Prmpt=[Prmpt,' (Default: ',num2str(DefItem),')'];
            end
        case 'mq'
            for i = 1:nLabels, 
                cprintf('*black','\t%2d | ', i);
                fprintf('%s\n',Labels(i,:)), 
            end
            cprintf('*black','\t<==| Back/Quit [Q]\n');
            Prmpt = ['Menu choice (1-',int2str(nLabels),'/Q)'];
            if DefItem
                Prmpt=[Prmpt,' (Default: ',num2str(DefItem),')'];
            else
                Prmpt=[Prmpt,' (Default: Q)'];
            end
            nLabels = nLabels+1;
            DefItem = nLabels;
            if iscell(Values)
                Values(end+1) = {'BACK'};
            else
                Values(end+1,:) = 0;
            end
    end
    fprintf('%s \n',repmat('_',1,LabelWidth+15));
    %-Ask for user response
    %-------------------------------------------------------
    if nLabels==1
        %-Only one choice - auto-pick & display
        k = 1;
        cprintf('red','Menu choice: 1 - %s\t(only option)',Labels)
    else
        switch lower(Type)
            case 'm'
               k = input(sprintf(' %s ? ',Prmpt));
            case 'mq'
               k = input(sprintf(' %s ? ',Prmpt),'s');
               if strcmpi(k,'q') || isempty(k)
                   k = numel(Values);
               else
                   k = str2double(k);
               end
        end
        if DefItem & isempty(k), k=DefItem; end
        while isempty(k) || ~any([1:nLabels]==k)
            if ~isempty(k),cprintf('_red','\n!Selected menu item does not exist \n'),end
            k = input([fprintf('%s',Prmpt) ' ? ']);
            if DefItem & isempty(k), k=DefItem; end
        end
    end
    fprintf('\n')

	p = Values(k,:); if ischar(p), p=deblank(p); end

otherwise, error('unrecognised type')
end % (switch lower(Type) within case {'b','b|','y/n'})

%-Log the transaction & return response
%-----------------------------------------------------------------------
if exist('spm_log')==2
	if iscellstr(Prompt), Prompt=Prompt{1}; end
	spm_log([mfilename,' : ',Prompt,': (',deblank(Labels(k,:)),')'],p); end
varargout = {p,YPos};

case {'m!','b!'}                          %-GUI PullDown/Buttons utility
%=======================================================================
% H = nk_input(Prompt,YPos,'p',Labels,cb,UD,XCB)
%-Condition arguments
if nargin<7, XCB    = 0;  else, XCB    = varargin{7}; end
if nargin<6, UD     = []; else, UD     = varargin{6}; end
if nargin<5, cb     = ''; else, cb     = varargin{5}; end
if nargin<4, Labels = []; else, Labels = varargin{4}; end

if CmdLine, error('Can''t do CmdLine GUI utilities!'), end
if isempty(cb), cb = 'disp(''(CallBack not set)'')'; end
if ischar(cb), cb = cellstr(cb); end
if length(cb)>1 & strcmp(lower(Type),'m!'), XCB=1; end

if iscellstr(Labels), Labels=char(Labels); end
%-Convert Labels "option" string to string matrix if required
if any(Labels=='|')
	OptStr=Labels;
	BarPos=find([OptStr=='|',1]);
	Labels=OptStr(1:BarPos(1)-1);
	for Bar = 2:sum(OptStr=='|')+1
		Labels=strvcat(Labels,OptStr(BarPos(Bar-1)+1:BarPos(Bar)-1));
	end
end

%-Check #CallBacks
if ~( length(cb)==1 | (length(cb)==size(Labels,1)) )
	error('Labels & Callbacks size mismatch'), end


%-Draw Prompt
%-----------------------------------------------------------------------
Tag = ['GUIinput_',int2str(YPos)];			%-Tag for widgets

if ~isempty(Prompt)
	uicontrol(Finter,'Style','Text',...
		'String',Prompt,...
		'Tag',Tag,...
		'HorizontalAlignment','Right',...
        'BackgroundColor',COLOUR,...
		'Position',PRec)
	Rec = RRec;
else
	Rec = QRec;
end


%-Sort out UserData for extended callbacks (handled by nk_input('!m_cb')
%-----------------------------------------------------------------------
if XCB, if iscell(UD), UD={UD}; end, UD = struct('UD',UD,'cb',{cb}); end


%-Draw PullDown or Buttons
%-----------------------------------------------------------------------
switch lower(Type), case 'm!'
	if XCB, UD.cb=cb; cb = {'nk_input(''!m_cb'')'}; end
	H = uicontrol(Finter,'Style','PopUp',...
		'HorizontalAlignment','Left',...
		'ForegroundColor','k',...
		'BackgroundColor',COLOUR,...
		'String',Labels,...
		'Tag',Tag,...
		'UserData',UD,...
		'CallBack',char(cb),...
		'Position',Rec);

case 'b!'
	nLabels = size(Labels,1);
	dX = Rec(3)/nLabels;

	H = [];
	for i=1:nLabels
		if length(cb)>1, tcb=cb(i); else, tcb=cb; end
		if XCB, UD.cb=tcb; tcb = {'nk_input(''!m_cb'')'}; end
		h = uicontrol(Finter,'Style','Pushbutton',...
			'String',deblank(Labels(i,:)),...
			'ToolTipString','',...
			'Tag',Tag,...
			'UserData',UD,...
			'BackgroundColor',COLOUR,...
			'Callback',char(tcb),...
			'Position',[Rec(1)+(i-1)*dX+1 ...
					Rec(2) dX-2 Rec(4)]);
		H = [H,h];
	end


end


%-Bring window to fore & jump pointer to menu widget
[PLoc,cF] = nk_input('!PointerJump',RRec,Finter);

varargout = {H};

case {'d','d!'}                                        %-Display message
%=======================================================================
%-Condition arguments
if nargin<4, Label=''; else, Label=varargin{4}; end

if CmdLine & strcmp(lower(Type),'d')
	fprintf('\n     +-%s%s+',Label,repmat('-',1,57-length(Label)))
	Prompt = [Prompt,' '];
	while length(Prompt)>0
		tmp = length(Prompt);
		if tmp>56, tmp=min([max(find(Prompt(1:56)==' ')),56]); end
		fprintf('\n     | %s%s |',Prompt(1:tmp),repmat(' ',1,56-tmp))
		Prompt(1:tmp)=[];
	end
	fprintf('\n     +-%s+\n',repmat('-',1,57))
%  elseif ~CmdLine
%  	if ~isempty(Label), Prompt = [Label,': ',Prompt]; end
%  	figure(Finter)
%  	%-Create text axes and edit control objects
%  	%---------------------------------------------------------------
%  	h = uicontrol(Finter,'Style','Text',...
%  		'String',Prompt(1:min(length(Prompt),56)),...
%  		'FontWeight','bold',...
%  		'Tag',['GUIinput_',int2str(YPos)],...
%  		'HorizontalAlignment','Left',...
%  		'ForegroundColor','k',...
%          'BackgroundColor',COLOUR,...
%  		'UserData',Prompt,...
%  		'Position',QRec);
%  	if length(Prompt)>56
%  		pause(1)
%  		set(h,'ToolTipString',Prompt)
%  		nk_input('!dScroll',h)
%  		uicontrol(Finter,'Style','PushButton','String','>',...
%  			'ToolTipString','press to scroll message',...
%  			'Tag',['GUIinput_',int2str(YPos)],...
%  			'UserData',h,...
%  			'CallBack',[...
%  			 'set(gcbo,''Visible'',''off''),',...
%  			 'nk_input(''!dScroll'',get(gcbo,''UserData'')),',...
%  			 'set(gcbo,''Visible'',''on'')'],...
%              'BackgroundColor',COLOUR,...
%  			'Position',[QRec(1)+QRec(3)-10,QRec(2),15,QRec(4)]);
%  	end
end
if nargout>0, varargout={[],YPos}; end



%=======================================================================
% U T I L I T Y   F U N C T I O N S 
%=======================================================================

case '!colour'
%=======================================================================
% colour = nk_input('!Colour')
varargout = {COLOUR};


case '!icond'
%=======================================================================
% [iCond,msg] = nk_input('!iCond',str,n,m)
% Parse condition indicator spec strings:
%	'2 3 2 3', '0 1 0 1', '2323', '0101', 'abab', 'R A R A'
if nargin<4, m=Inf; else, m=varargin{4}; end
if nargin<3, n=NaN; else, n=varargin{3}; end
if any(isnan(n(:)))
	n=Inf;
elseif (length(n(:))==2 & ~any(n==1)) | length(n(:))>2
	error('condition input can only do vectors')
end
if nargin<2, i=''; else, i=varargin{2}; end
if isempty(i), varargout={[],'empty input'}; return, end
msg = ''; i=i(:)';

if ischar(i)
	if i(1)=='0' & all(ismember(unique(i(:)),setstr(abs('0'):abs('9'))))
		%-Leading zeros in a digit list
		msg = sprintf('%s expanded',i);
		z = min(find([diff(i=='0'),1]));
		i = [zeros(1,z), nk_input('!iCond',i(z+1:end))'];
	else
		%-Try an eval, for functions & string #s
		i = evalin('base',['[',i,']'],'i');
	end
end

if ischar(i)
	%-Evaluation error from above: see if it's an 'abab' or 'a b a b' type:
	[c,null,i] = unique(lower(i(~isspace(i))));
	if all(ismember(c,setstr(abs('a'):abs('z'))))
		%-Map characters a-z to 1-26, but let 'r' be zero (rest)
		tmp = c-'a'+1; tmp(tmp=='r'-'a'+1)=0;
		i   = tmp(i);
		msg = [sprintf('[%s] mapped to [',c),...
			sprintf('%d,',tmp(1:end-1)),...
			sprintf('%d',tmp(end)),']'];
	else
		i = '!'; msg = 'evaluation error';
	end
elseif ~all(floor(i(:))==i(:))
	i = '!'; msg = 'must be integers';
elseif length(i)==1 & prod(n)>1
	msg = sprintf('%d expanded',i);
	i = floor(i./10.^[floor(log10(i)+eps):-1:0]);
	i = i-[0,10*i(1:end-1)];
end

%-Check size of i & #conditions
if ~ischar(i), [i,msg] = sf_SzChk(i,n,msg); end
if ~ischar(i) & isfinite(m) & length(unique(i))~=m
	i = '!'; msg = sprintf('%d conditions required',m);
end

varargout = {i,msg};


case '!inptconmen'
%=======================================================================
% hM = nk_input('!InptConMen',Finter,H)
if nargin<3, H=[]; else, H=varargin{3}; end
if nargin<2, varargout={[]}; else, Finter=varargin{2}; end
hM = uicontextmenu('Parent',Finter);
uimenu(hM,'Label','help on nk_input',...
	'CallBack','spm_help(''nk_input.m'')')
if ConCrash
	uimenu(hM,'Label','crash out','Separator','on',...
		'CallBack','delete(get(gcbo,''UserData''))',...
		'UserData',[hM,H])
end

set(Finter,'UIContextMenu',hM)

varargout={hM};


case '!cmdline'
%=======================================================================
% [CmdLine,YPos] = nk_input('!CmdLine',YPos)
%-Sorts out whether to use CmdLine or not & canonicalises YPos
if nargin<2, YPos=''; else, YPos=varargin{2}; end
if isempty(YPos), YPos='+1'; end

CmdLine = [];

%-Special YPos specifications
if ischar(YPos)
	if(YPos(1)=='!'), CmdLine=0; YPos(1)=[]; end
elseif YPos==0
	CmdLine=1;
elseif YPos<0
	CmdLine=0;
	YPos=-YPos;
end

CmdLine = spm('CmdLine',CmdLine);
if CmdLine, YPos=0; end

varargout = {CmdLine,YPos};


case '!getwin'
%=======================================================================
% Finter = nk_input('!GetWin',F)
%-Locate (or create) figure to work in (Don't use 'Tag'ged figs)
if nargin<2, F='Interactive'; else, F=varargin{2}; end
Finter = spm_figure('FindWin',F);
if isempty(Finter)
	if any(get(0,'Children'))
		if isempty(get(gcf,'Tag')), Finter = gcf;
		else, Finter = spm('CreateIntWin'); end
	else, Finter = spm('CreateIntWin'); end
end
varargout = {Finter};


case '!pointerjump'
%=======================================================================
% [PLoc,cF] = nk_input('!PointerJump',RRec,F,XDisp)
%-Raise window & jump pointer over question
if nargin<4, XDisp=[]; else, XDisp=varargin{4}; end
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, error('Insufficient arguments'), else, RRec=varargin{2}; end
F = spm_figure('FindWin',F);
PLoc = get(0,'PointerLocation');
cF   = get(0,'CurrentFigure');
if ~isempty(F)
	figure(F)
	FRec = get(F,'Position');
	if isempty(XDisp), XDisp=RRec(3)*4/5; end
	if PJump, set(0,'PointerLocation',...
		floor([(FRec(1)+RRec(1)+XDisp), (FRec(2)+RRec(2)+RRec(4)/3)]));
	end
end
varargout = {PLoc,cF};


case '!pointerjumpback'
%=======================================================================
% nk_input('!PointerJumpBack',PLoc,cF)
%-Replace pointer and reset CurrentFigure back
if nargin<4, cF=[]; else, F=varargin{3}; end
if nargin<2, error('Insufficient arguments'), else, PLoc=varargin{2}; end
if PJump, set(0,'PointerLocation',PLoc), end
cF = spm_figure('FindWin',cF);
if ~isempty(cF), set(0,'CurrentFigure',cF); end


case '!prntprmpt'
%=======================================================================
% nk_input('!PrntPrmpt',Prompt,TipStr,Title)
%-Print prompt for CmdLine questioning
if nargin<5, lbw  = 72; else, lbw = varargin{5}+15; end
if nargin<4, Title  = ''; else, Title  = varargin{4}; end
if nargin<3, TipStr = ''; else, TipStr = varargin{3}; end
if nargin<2, Prompt = ''; else, Prompt = varargin{2}; end
if isempty(Prompt), Prompt='Enter an expression'; end

Prompt = cellstr(Prompt);

if ~isempty(TipStr)
  tmp    = 8 + length(Prompt{end}) + length(TipStr);
  if tmp < 62
    TipStr = sprintf('%s(%s)',repmat(' ',1,lbw-tmp),TipStr);
  else
    TipStr = sprintf('\n%s(%s)',repmat(' ',1,max(0,lbw-length(TipStr))),TipStr);
  end
end

if isempty(Title)
	fprintf('\n%s\n',repmat('~',1,lbw))
else
	fprintf('\n= ');cprintf('*black',' %s \n',Title)
    fprintf('%s',repmat('~',1,lbw-length(Title)-3));
end
cprintf(sprintf('*[%g,%g,%g]',NMinfo.clmenu(1),NMinfo.clmenu(2),NMinfo.clmenu(3)),'\t%s ',Prompt{1})
for i=2:prod(size(Prompt)), fprintf('\n\t%s',Prompt{i}), end
fprintf('%s\n%s\n',TipStr,repmat('~',1,lbw))

case '!inputrects'
%=======================================================================
% [Frec,QRec,PRec,RRec,Sz,Se] = nk_input('!InputRects',YPos,rec,F)
if nargin<4, F='Interactive'; else, F=varargin{4}; end
if nargin<3, rec=''; else, rec=varargin{3}; end
if nargin<2, YPos=1; else, YPos=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('Figure not found'), end

Units = get(F,'Units');
set(F,'Units','pixels')
FRec = get(F,'Position');
set(F,'Units',Units);
Xdim = FRec(3); Ydim = FRec(4);

WS   = spm('WinScale');
Sz   = round(22*min(WS));	%-Height
Pd   = Sz/2;			%-Pad
Se   = 2*round(25*min(WS)/2);	%-Seperation
Yo   = round(2*min(WS));	%-Y offset for responses

a = 5.5/10;
y = Ydim - Se*YPos;
QRec   = [Pd            y         Xdim-2*Pd        Sz]; %-Question
PRec   = [Pd            y     floor(a*Xdim)-2*Pd   Sz]; %-Prompt
RRec   = [ceil(a*Xdim)  y+Yo  floor((1-a)*Xdim)-Pd Sz]; %-Response
% MRec = [010           y         Xdim-50          Sz]; %-Menu PullDown
% BRec = MRec + [Xdim-50+1, 0+1, 50-Xdim+30, 0];        %-Menu PullDown OK butt

if ~isempty(rec)
	varargout = {eval(rec)};
else
	varargout = {FRec,QRec,PRec,RRec,Sz,Se};
end


case '!deleteinputobj'
%=======================================================================
% nk_input('!DeleteInputObj',F)
if nargin<2, F='Interactive'; else, F=varargin{2}; end
h = nk_input('!FindInputObj',F);
delete(h(h>0))


case {'!currentpos','!findinputobj'}
%=======================================================================
% [CPos,hCPos] = nk_input('!CurrentPos',F)
% h            = nk_input('!FindInputObj',F)
% hPos contains handles: Columns contain handles corresponding to Pos
if nargin<2, F='Interactive'; else, F=varargin{2}; end
F = spm_figure('FindWin',F);

%-Find tags and YPos positions of 'GUIinput_' 'Tag'ged objects
H    = [];
YPos = [];
for h = get(F,'Children')'
	tmp = get(h,'Tag');
	if ~isempty(tmp)
		if strcmp(tmp(1:min(length(tmp),9)),'GUIinput_')
			H    = [H, h];
			YPos = [YPos, eval(tmp(10:end))];
		end
	end
end

switch lower(Type), case '!findinputobj'
	varargout = {H};
case '!currentpos'
	if nargout<2
		varargout = {max(YPos),[]};
	elseif isempty(H)
		varargout = {[],[]};
	else
		%-Sort out 
		tmp     = sort(YPos);
		CPos    = tmp(find([1,diff(tmp)]));
		nPos    = length(CPos);
		nPerPos = diff(find([1,diff(tmp),1]));
		hCPos   = zeros(max(nPerPos),nPos);
		for i = 1:nPos
			hCPos(1:nPerPos(i),i) = H(YPos==CPos(i))';
		end
		varargout = {CPos,hCPos};
	end
end

case '!nextpos'
%=======================================================================
% [NPos,CPos,hCPos] = nk_input('!NextPos',YPos,F,CmdLine)
%-Return next position to use
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, YPos='+1'; else, YPos=varargin{2}; end
if nargin<4, [CmdLine,YPos]=nk_input('!CmdLine',YPos);
	else, CmdLine=varargin{4}; end

F = spm_figure('FindWin',F);

%-Get current positions
if nargout<3
	CPos = nk_input('!CurrentPos',F);
	hCPos = [];
else
	[CPos,hCPos] = nk_input('!CurrentPos',F);
end

if CmdLine
	NPos = 0;
else
	MPos = nk_input('!MaxPos',F);
	if ischar(YPos)
		%-Relative YPos
		%-Strip any '!' prefix from YPos
		if(YPos(1)=='!'), YPos(1)=[]; end
		if strncmp(YPos,'_',1)
			%-YPos='_' means bottom
			YPos=eval(['MPos+',YPos(2:end)],'MPos');
		else
			YPos = max([0,CPos])+eval(YPos);
		end
	else
		%-Absolute YPos
		YPos=abs(YPos);
	end
	NPos = min(max(1,YPos),MPos);
end
varargout = {NPos,CPos,hCPos};

case '!setnextpos'
%=======================================================================
% NPos = nk_input('!SetNextPos',YPos,F,CmdLine)
%-Set next position to use
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, YPos='+1'; else, YPos=varargin{2}; end
if nargin<4, [CmdLine,YPos]=nk_input('!CmdLine',YPos);
	else CmdLine=varargin{4}; end

%-Find out which Y-position to use
[NPos,CPos,hCPos] = nk_input('!NextPos',YPos,F,CmdLine);

%-Delete any previous inputs using positions NPos and after
if any(CPos>=NPos), h=hCPos(:,CPos>=NPos); delete(h(h>0)), end

varargout = {NPos};

case '!maxpos'
%=======================================================================
% MPos = nk_input('!MaxPos',F,FRec3)
%
if nargin<3
	if nargin<2, F='Interactive'; else, F=varargin{2}; end
	F = spm_figure('FindWin',F);
	if isempty(F)
		FRec3=spm('WinSize','Interactive')*[0;0;0;1];
	else
		%-Get figure size
		Units = get(F,'Units');
		set(F,'Units','pixels')
		FRec3 = get(F,'Position')*[0;0;0;1];
		set(F,'Units',Units);
	end
end

Se   = round(25*min(spm('WinScale')));
MPos = floor((FRec3-5)/Se);

varargout = {MPos};


case '!editablekeypressfcn'
%=======================================================================
% nk_input('!EditableKeyPressFcn',h,ch,hPrmpt)
if nargin<2, error('Insufficient arguments'), else, h=varargin{2}; end
if isempty(h), set(gcbf,'KeyPressFcn','','UserData',[]), return, end
if nargin<3, ch=get(get(h,'Parent'),'CurrentCharacter'); else, ch=varargin{3};end
if nargin<4, hPrmpt=get(h,'UserData'); else, hPrmpt=varargin{4}; end

tmp = get(h,'String');
if isempty(tmp), tmp=''; end
if iscellstr(tmp) & length(tmp)==1; tmp=tmp{:}; end

if isempty(ch)					%- shift / control / &c. pressed
	return
elseif any(abs(ch)==[32:126])			%-Character
	if iscellstr(tmp), return, end
	tmp = [tmp, ch];
elseif abs(ch)==21				%- ^U - kill
	tmp = '';
elseif any(abs(ch)==[8,127])			%-BackSpace or Delete
	if iscellstr(tmp), return, end
	if length(tmp), tmp(length(tmp))=''; end
elseif abs(ch)==13				%-Return pressed
	if ~isempty(tmp)
		set(hPrmpt,'UserData',get(h,'String'))
	end
	return
else
	%-Illegal character
	return
end
set(h,'String',tmp)


case '!buttonkeypressfcn'
%=======================================================================
% nk_input('!ButtonKeyPressFcn',h,Keys,DefItem,ch)
%-Callback for KeyPress, to store valid button # in UserData of Prompt,
% DefItem if (DefItem~=0) & return (ASCII-13) is pressed

%-Condition arguments
if nargin<2, error('Insufficient arguments'), else, h=varargin{2}; end
if isempty(h), set(gcf,'KeyPressFcn','','UserData',[]), return, end
if nargin<3, error('Insufficient arguments'); else, Keys=varargin{3}; end
if nargin<4, DefItem=0; else, DefItem=varargin{4}; end
if nargin<5, ch=get(gcf,'CurrentCharacter'); else, ch=varargin{5}; end

if isempty(ch)
	%- shift / control / &c. pressed
	return
elseif (DefItem & ch==13)
	But = DefItem;
else
	But = find(lower(ch)==lower(Keys));
end
if ~isempty(But), set(h,'UserData',But), end


case '!pulldownkeypressfcn'
%=======================================================================
% nk_input('!PullDownKeyPressFcn',h,ch,DefItem)
if nargin<2, error('Insufficient arguments'), else, h=varargin{2}; end
if isempty(h), set(gcf,'KeyPressFcn',''), return, end
if nargin<3, ch=get(get(h,'Parent'),'CurrentCharacter'); else, ch=varargin{3};end
if nargin<4, DefItem=get(h,'UserData'); else, ch=varargin{4}; end

Pmax = get(h,'Max');
Pval = get(h,'Value');

if Pmax==1, return, end

if isempty(ch)
	%- shift / control / &c. pressed
	return
elseif abs(ch)==13
	if Pval==1
		if DefItem, set(h,'Value',max(2,min(DefItem+1,Pmax))), end
	else
		set(h,'UserData','Selected')
	end
elseif any(ch=='bpu')
	%-Move "b"ack "u"p to "p"revious entry
	set(h,'Value',max(2,Pval-1))
elseif any(ch=='fnd')
	%-Move "f"orward "d"own to "n"ext entry
	set(h,'Value',min(Pval+1,Pmax))
elseif any(ch=='123456789')
	%-Move to entry n
	set(h,'Value',max(2,min(eval(ch)+1,Pmax)))
else
	%-Illegal character
end


case '!m_cb'     %-CallBack handler for extended CallBack 'p'ullDown type
%=======================================================================
% nk_input('!m_cb')

%-Get PopUp handle and value
h   = gcbo;
n   = get(h,'Value');

%-Get PopUp's UserData, check cb and UD fields exist, extract cb & UD
tmp = get(h,'UserData');
if ~(isfield(tmp,'cb') & isfield(tmp,'UD'))
	error('Invalid UserData structure for nk_input extended callback')
end
cb  = tmp.cb;
UD  = tmp.UD;

%-Evaluate appropriate CallBack string (ignoring any return arguments)
% NB: Using varargout={eval(cb{n})}; gives an error if the CallBack 
% has no return arguments!
if length(cb)==1, eval(char(cb)); else, eval(cb{n}); end


case '!dscroll'
%=======================================================================
% nk_input('!dScroll',h,Prompt)
%-Scroll text in object h
if nargin<2, return, else, h=varargin{2}; end
if nargin<3, Prompt = get(h,'UserData'); else, Prompt=varargin{3}; end
tmp = Prompt;
if length(Prompt)>56
	while length(tmp)>56
		tic, while(toc<0.1), pause(0.05), end
		tmp(1)=[];
		set(h,'String',tmp(1:min(length(tmp),56)))
	end
	pause(1)
	set(h,'String',Prompt(1:min(length(Prompt),56)))
end


otherwise
%=======================================================================
error(['Invalid type/action: ',Type])


%=======================================================================
end % (case lower(Type))


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function [Keys,Labs] = sf_labkeys(Labels)
%=======================================================================
%-Make unique character keys for the Labels, ignoring case
if nargin<1, error('insufficient arguments'), end
if iscellstr(Labels), Labels = char(Labels); end
if isempty(Labels), Keys=''; Labs=''; return, end

Keys=Labels(:,1)';
nLabels = size(Labels,1);
if any(~diff(abs(sort(lower(Keys)))))
	if nLabels<10
		Keys = sprintf('%d',[1:nLabels]);
	elseif nLabels<=26
		Keys = sprintf('%c',abs('a')+[0:nLabels-1]);
	else
		error('Too many buttons!')
	end
	Labs = Labels;
else
	Labs = Labels(:,2:end);
end


function [p,msg] = sf_eEval(str,Type,n,m,l)
%=======================================================================
%-Evaluation and error trapping of typed input
if nargin<5, l=[]; end
if nargin<4, m=[]; end
if nargin<3, n=[]; end
if nargin<2, Type='e'; end
if nargin<1, str=''; end
if isempty(str) 
    if strcmp(Type,'sq')
        p=[]; return
    else
        p='!'; msg='empty input'; return, 
    end
end
switch lower(Type)
case {'s','sq'}
	p = str; msg = '';
case 'e'
	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	else
		[p,msg] = sf_SzChk(p,n);
	end
case 'n'
	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	elseif any(floor(p(:))~=p(:)|p(:)<1)|~isreal(p)
		p='!'; msg='natural number(s) required';
	elseif ~isempty(m) & any(p(:)>m)
		p='!'; msg=['max value is ',num2str(m)];
	else
		[p,msg] = sf_SzChk(p,n);
	end
case 'w'
	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	elseif any(floor(p(:))~=p(:)|p(:)<0)|~isreal(p)
		p='!'; msg='whole number(s) required';
	elseif ~isempty(m) & any(p(:)>m)
		p='!'; msg=['max value is ',num2str(m)];
    elseif ~isempty(l) & any(p(:)<l)
        p='!'; msg=['min value is ',num2str(l)];
	else
		[p,msg] = sf_SzChk(p,n);
	end
case 'i'
	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	elseif any(floor(p(:))~=p(:))|~isreal(p)
		p='!'; msg='integer(s) required';
	else
		[p,msg] = sf_SzChk(p,n);
	end
case 'p'
	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	elseif length(setxor(p(:)',m))
		p='!'; msg='invalid permutation';
	else
		[p,msg] = sf_SzChk(p,n);
	end
case 'r'
	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	elseif ~isreal(p)
		p='!'; msg='real number(s) required';
	elseif ~isempty(m) & ( max(p)>max(m) | min(p)<min(m) )
		p='!'; msg=sprintf('real(s) in [%g,%g] required',min(m),max(m));
	else
		[p,msg] = sf_SzChk(p,n);
	end
case 'c'
	if isempty(m), m=Inf; end
	[p,msg] = nk_input('!iCond',str,n,m);
case 'x'
	X = m;			%-Design matrix/space-structure
	if isempty(n), n=1; end

	%-Sort out contrast matrix dimensions (contrast vectors in rows)
	if length(n)==1, n=[n,Inf]; else, n=reshape(n(1:2),1,2); end
	if ~isempty(X)		% - override n(2) w/ design column dimension
		n(2) = spm_SpUtil('size',X,2);
	end

	p = evalin('base',['[',str,']'],'''!''');
	if ischar(p)
		msg = 'evaluation error';
	else
		if isfinite(n(2)) & size(p,2)<n(2)
			tmp = n(2) -size(p,2);
			p   = [p, zeros(size(p,1),tmp)];
			if size(p,1)>1, str=' columns'; else, str='s'; end
			msg = sprintf('right padded with %d zero%s',tmp,str);
		else
			msg = '';
		end

		if size(p,2)>n(2)
			p='!'; msg=sprintf('too long - only %d prams',n(2));
		elseif isfinite(n(1)) & size(p,1)~=n(1)
			p='!';
			if n(1)==1, msg='vector required';
			    else, msg=sprintf('%d contrasts required',n(1)); end
		elseif ~isempty(X) & ~spm_SpUtil('allCon',X,p')
			p='!'; msg='invalid contrast';
		end
	end

otherwise
	error('unrecognised type');
end


function str = sf_SzStr(n,l)
%=======================================================================
%-Size info string constuction
if nargin<2, l=0; else, l=1; end
if nargin<1, error('insufficient arguments'), end
if isempty(n), n=NaN; end
n=n(:); if length(n)==1, n=[n,1]; end, dn=length(n);
if any(isnan(n)) | (prod(n)==1 & dn<=2) | (dn==2 & min(n)==1 & isinf(max(n)))
	str = ''; lstr = '';
elseif dn==2 & min(n)==1
	str = sprintf('[%d]',max(n));	lstr = [str,'-vector'];
elseif dn==2 & sum(isinf(n))==1
	str = sprintf('[%d]',min(n));	lstr = [str,'-vector(s)'];
else
	str=''; for i = 1:dn
		if isfinite(n(i)), str = sprintf('%s,%d',str,n(i));
		else, str = sprintf('%s,*',str); end
	end
	str = ['[',str(2:end),']'];	lstr = [str,'-matrix'];
end
if l, str=sprintf('\t%s',lstr); else, str=[str,' ']; end


function [p,msg] = sf_SzChk(p,n,msg)
%=======================================================================
%-Size checking
if nargin<3, msg=''; end
if nargin<2, n=[]; end, if isempty(n), n=NaN; else, n=n(:)'; end
if nargin<1, error('insufficient arguments'), end

if ischar(p) | any(isnan(n(:))), return, end
if length(n)==1, n=[n,1]; end

dn = length(n);
sp = size(p);
dp = ndims(p);

if dn==2 & min(n)==1
	%-[1,1], [1,n], [n,1], [1,Inf], [Inf,1] - vector - allow transpose
	%---------------------------------------------------------------
	i = min(find(n==max(n)));
	if n(i)==1 & max(sp)>1
		p='!'; msg='scalar required';
	elseif ndims(p)~=2 | ~any(sp==1) | ( isfinite(n(i)) & max(sp)~=n(i) )
		%-error: Not2D | not vector | not right length
		if isfinite(n(i)), str=sprintf('%d-',n(i)); else, str=''; end
		p='!'; msg=[str,'vector required'];
	elseif sp(i)==1 & n(i)~=1
		p=p'; msg=[msg,' (input transposed)'];
	end

elseif dn==2 & sum(isinf(n))==1
	%-[n,Inf], [Inf,n] - n vector(s) required - allow transposing
	%---------------------------------------------------------------
	i = find(isfinite(n));
	if ndims(p)~=2 | ~any(sp==n(i))
		p='!'; msg=sprintf('%d-vector(s) required',min(n));
	elseif sp(i)~=n
		p=p'; msg=[msg,' (input transposed)'];
	end	

else
	%-multi-dimensional matrix required - check dimensions
	%---------------------------------------------------------------
	if ndims(p)~=dn | ~all( size(p)==n | isinf(n) )
		p = '!'; msg='';
		for i = 1:dn
			if isfinite(n(i)), msg = sprintf('%s,%d',msg,n(i));
			else, msg = sprintf('%s,*',msg); end
		end
		msg = ['[',msg(2:end),']-matrix required'];
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        UIFOCUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function uifocus(h)
% MATLAB R14 NO LONGER KEEPS FOCUS ON THE UICONTROL, UNLESS IT IS
% EXPLICITLY SPECIFIED, OTHERWISE CONTROL IS GIVEN TO THE FIGURE'S
% KEYPRESS FUNCTION. THEREFORE WE MUST EXPLICITLY SET THE FOCUS TO
% A UICONTROL OR THE USER IS FORCED TO CLICK INTO THE CONTROL TO
% ENTER DATA OR PRESS A BUTTON. MATLAB 7 HAS EXTENDED THE FUNCTIONALITY
% OF THE UICONTROL FUNCTION. UICONTROL(HANDLE) SETS THE FOCUS
% TO THE DESIGNATED UIONTROL'S HANDLE.
% -DRG 05/1/13
%------------------------------------------------------------------
a = ver('Matlab');
if str2num(a.Version) > 7
    uicontrol(h);
end
return
%-------------------------------------------------------------------

