function out = tabselectionfcn(hFig,tag,count,action)
% TABSELECTIONFCN  allows to switch and enable/disable the tabpanels programatically
% Using the TABSELECTIONFCN it is possible to switch, enable and disable
% the tabs (that were created using the TabPanel Constructor) programatically.
%
% usage:
% 1. Select Tabpanel
% ------------------
% TABSELECTIONFCN(<hFig>,<TabTag>,<tabnumber>)
%     <hFig>      the handle(!) of the Figure (which contains the tabpanel)
%                 and not the name of the figure file.
%     <TabTag>    the Tag name of the tabpanel
%     <tabnumber> The number of the tabpanel or the tab string
%
% 2. Enable/Disable tabpanel(s)
% -----------------------------
% TABSELECTIONFCN(<hFig>,<TabTag>,<tabnumber(s)>,'on/off')
%     <tabnumber(s)> can be a scalar or vector of indices.
%
% Note: You can not disable the tabpanel that is currently active.
%
% 3. Obtain current state of all panels
% -------------------------------------
% state = TABSELECTIONFCN(<hFig>,<TabTag>)
%     where "state" is a vector contains:
%      1 - active tab panel
%      0 - tab is enable, but not active
%     -1 - tab is disabled
%
% See also TABPANEL.
%
%   Version: v1.3
%      Date: 2010/06/18 12:00:00
%   (c) 2008 By Elmar Tarajan [MCommander@gmx.de]
%
%
%  <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/4137">my other tools from File Exchange</a>
%
%  <a href="https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=Z8NDTASMAYB8W">PayPal Donate</a> :-)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% I LOVE MATLAB! You too? :) %%%