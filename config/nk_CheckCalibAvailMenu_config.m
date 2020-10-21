function [menustr, menuact] = nk_CheckCalibAvailMenu_config(menustr, menuact, calibuse)
global CALIBAVAIL

if CALIBAVAIL
    
    if exist('calibuse','var') && ~isempty(calibuse) && calibuse 
        calibstr = 'yes';
    else
        calibstr = 'no';
    end
    
    menustr = [sprintf('Use calibration data [ %s ]|',calibstr) menustr]; 
    menuact = [ 1000 menuact ];
    
end

