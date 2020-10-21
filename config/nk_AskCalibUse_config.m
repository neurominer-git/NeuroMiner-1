function CALIBUSE = nk_AskCalibUse_config(ASKSTR, CALIBUSE)

CALIBUSE = nk_input(sprintf('Do you want to use the calibration sample for %s',ASKSTR), 0, 'yes|no', [1,2], CALIBUSE);
