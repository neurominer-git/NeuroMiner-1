function [IO, groupmode_str] = SetFileFilter(IO,groupmode,datasource)
global SPMAVAIL

switch groupmode
    case -1
        groupmode_str = 'SPM.mat file';
        IO.filt = 'SPM.mat$'; IO.M_filestr = 'SPM.mat';
        IO.spacedef_filt = '.*\.nii$|.*\.img$';
    case 0
        groupmode_str = 'image files';
        switch datasource 
            case 'nifti'
                IO.filt = '.*\.nii$|.*\.img$';
                IO.spacedef_filt = '.*\.nii$|.*\.img$';
            case 'surf'
                IO.filt = '.*\.mgh$|.*\.mgz$|.*\.gii$';
                IO.spacedef_filt = '.*\.mgh$|.*\.mgz$|.*\.gii$';
            case 'vector'
                IO.filt = '.*\.csv$|.*\.txt$|.*\.xml$';
                IO.spacedef_filt = [];
        end
    case 1
        groupmode_str = 'the MATLAB workspace';
        IO.filt = []; IO.M_filestr = 'MATLAB';
        IO.spacedef_filt = [];
    case 2
        groupmode_str = 'the MATLAB file';
        if SPMAVAIL
            IO.filt = '.*\.mat$'; 
        else
            IO.filt = {'*.mat'}; 
        end
        IO.M_filestr = 'MATLAB';
        IO.spacedef_filt = [];
    case 3
        groupmode_str = 'the text file';
        if SPMAVAIL
            IO.filt = '.*\.txt$|.*\.csv';
        else
            IO.filt = {'*.txt;*.csv'};
        end
         IO.M_filestr = 'CSV/TXT';
        IO.spacedef_filt = [];
    case 4
        groupmode_str = 'the spreadsheet';
        if SPMAVAIL
            IO.filt = '.*\.xlsx$|.*\.xls$|.*\.ods'; 
        else
            IO.filt = {'*.xlsx;*.xls;*.ods'}; 
        end
        IO.M_filestr = 'spreadsheet';
        IO.spacedef_filt = [];
    case 5
        groupmode_str = 'the XML files';
        if SPMAVAIL
            IO.filt = '.*\.xml$'; 
        else
            IO.filt = {'*.xml'};
        end
            IO.M_filestr = 'XML';
        IO.spacedef_filt = [];
        
end