function tabpanel(figname,tag,action)
%TABPANEL  "TabPanel Constructor" offers the easiest way for creating tabpanels in MATLAB
%  Initialize Tabpanel using the the size and location:
%  >> tabpanel('filename.fig','tabpaneltag',[left bottom width height])
%     'tabpaneltag' define the tag name of the tabpanel. The rectangle defines 
%     the size and location (in pixel) of the tabpanel within the parent 
%     figure window.
%
%  Alternative (as usual):
%  1. Open the figure (where the tab panel needs to be integrated) in
%     GUIDE and create a Text-Object which reserves the place of the future
%     tabpanel.
%  2. Specify and note the unique tag-name (tabpaneltag) of the created 
%     text-object and the figure filename.
%  3. Start "TabPanel Constructor" as follows:
%        >> tabpanel('filename.fig','tabpaneltag')
%
%  Options:
%     a. activate "TabPanel Constructor" to edit an existing tabpanel:
%        >> tabpanel('filename.fig','tabpaneltag') 
%
%     b. remove tabpanel from GUI:
%        >> tabpanel('filename.fig','tabpaneltag','delete')
%
%     c. Querying the current state of tab panels.
%        The are two ways to obtain this info:
%        1) using the TABSELECTIONFCN
%        >> state = TABSELECTIONFCN(<hFig>,<TabTag>)
%            <hFig> is the handle(!) of the Figure (which contains the tabpanel)
%            and not the name of the figure file.
%
%        2) without using the TABSELECTIONFCN
%        >> state = get(<tabpanel-handle>,'value')
%            <tabpanel-handle> is the handle of the UI-Object, which was use
%            to initialize the tabpanel in the first step.
%
%        If you want to use this feature in tabpanel which was created using
%        the older TPC-Version so please follow the steps to update mentioned below.
%
%
%  !!! IMPORTANT NOTES !!!
%  1) While TPC is active, you can not use your figure as usual!
%     You should close TPC and restart your GUI application!
%
%  2) Follow the steps to update tabpanel to a current TPC version
%     1. Backup you FIG-File (and appropriate MATLAB-File)
%     2. Open the tabpanel with the current version
%     3. Add any new tab panel
%     4. Remove last created tab panel
%     5. Save and Close TPC
%
%  3) In case when the figure name has been changed you should just 
%     open the TPC and close it again.
%
%
% See also TABSELECTIONFCN.
%
%   Version: v2.8.1
%      Date: 2010/06/21 12:00:00
%   (c) 2005-2010 By Elmar Tarajan [MCommander@gmx.de]
%
%
%  <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/4137">my other tools from File Exchange</a>
%
%  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=Z8NDTASMAYB8W">PayPal Donate</a> :-)

%  History:
%   Version: v2.8.1
%   2010/06/25 fixed: error while adding "SelectionChangeFcn"
%  
%   Version: v2.8
%   2010/06/16 added: querying of the tab panel state
%   2010/06/10 added: adding already existing figures as tab
%   2010/06/09 added: "Close All" button while closing tabpanel
%   2010/05/19 fixed: error while removing tabpanel from figure
%   2010/05/09 some code improvements
%
%   Version: v2.7
%   2010/03/04 added: creating tabpanel using location and size
%   2010/03/04 added: adding already existing figures as tab (beta)
%   2010/01/15 fixed: Problems using TPC after renaming the main figure
%
%   Version: v2.6.4
%  *2010/03/04 todo: creating tabpanel using location and size
%   2009/09/16 improved: code for adding the "TabSelectionChange_Callback"
%
%   Version: v2.6.3
%  *2009/09/16 todo:  Allow to edit TabPanel-Tag after creating
%   2009/09/16 improved: tab-callback code
%  *2009/09/16 todo:  Improve code for adding the "TabSelectionChange_Callback"
%   2009/09/16 fixed: Error while saving figure which name has been changed
%   2009/08/12 added: adding already available figure as new panel (beta)
%   2009/07/24 fixed: Error while closing TPC in R2009a/b
%   2009/07/24 fixed: GUI-Restart while tab switching
%   2009/07/24 some code improvements
%
%   Version: v2.6.1
%   2008/09/10 added: "tabselectionfcn" for programatically tab switching
%   2008/08/03 added: TabSelectionChange_Callback
%   2008/07/02 supporting of nested tabpanels
%   2008/06/23 works now with R14/SP1/SP2/SP3 versions
%   2008/05/08 many code improvements
%   2008/04/16 improved look - using the mouse on "settings"-button
%   2008/03/28 some code improvements
%   2008/03/17 improved look of tabs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I LOVE MATLAB! You too? :) %%%