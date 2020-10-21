function display_writetable(tblstruct, filename, tblnames, sheetnames)
% =========================================================================
% FORMAT display_writetable(tblstruct, filename, tblnames, sheetnames)
% =========================================================================
% helper function of statistical comparison module in NM results manager
% Decides how to save results files to disk depending on the availability
% of Excel on OS. If Excel is not available the results files are saved as
% independent CSV files.
% =========================================================================
% (c) Nikolaos Koutsouleris, 08/2020
if ispc || ismac
    try
        objExcel = actxserver('Excel.Application');
        excelinstalled = true;

    catch
        fprintf('\nExcel not installed')
        excelinstalled=false;
    end
else
    fprintf('\nExcel not available on Linux OS')
    excelinstalled=false;
end

if excelinstalled
    % save results as excel spreadsheet
    fprintf('\nExcel installed on system,');
    for i=1:numel(sheetnames)
        fprintf('\nWriting sheet ''%s'' to excel file %s',sheetnames{i},filename); 
        writetable(tblstruct.(tblnames{i}),filename,'Sheet',sheetnames{i},'WriteRowNames',true);
    end
    try
        objExcel.Workbooks.Open(filename); 
        % MRZ: delete default sheets in freshly created .xls
        [~, sheets] = xlsfinfo(filename);
        sheetNames2remove = setdiff(sheets,sheetnames); 
        for i = 1:numel(sheetNames2remove)
            objExcel.ActiveWorkbook.Worksheets.Item(sheetNames2remove{i}).Delete;
        end
        objExcel.ActiveWorkbook.Save
        objExcel.ActiveWorkbook.Close(filename);
        delete(objExcel);
    end
else
    %... or as csv files
    [pth, nam, ext] = fileparts(filename); 
    for i=1:numel(sheetnames)
        filenames = fullfile(pth,[nam '_' sheetnames{i} ext]);
        fprintf('\nWriting %s to disk.',filenames)
        writetable(tblstruct.(tblnames{i}),filenames,'WriteRowNames',true);
    end 
  
end
fprintf('\nDone.\n')
