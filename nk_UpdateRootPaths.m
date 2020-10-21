function NM = nk_UpdateRootPaths(NM, AnalVec, NewRootDir)

if ~isfield(NM,'analysis')
   fprintf('\nNo analyses found in NM structure')
   return
end

if ~exist('NewRootDir','var') || isempty(NewRootDir) || ~exist(NewRootDir,'dir')
    NewRootDir = nk_DirSelector('Update analyses'' root paths');
end

if ~exist('AnalVec','var') || isempty(AnalVec)
    AnalVec = 1:numel(NM.analysis);
end

nA = numel(AnalVec);
failed = []; succeeded = []; cntf=1;cnts=1;
for j=1:nA
     try
         analind = AnalVec(j);
         if exist(NewRootDir,'dir')
            [~,analdir] = fileparts(NM.analysis{analind}.rootdir);
            NewAnalDir = fullfile(NewRootDir,analdir);
            % Adjust root directories
            NM.analysis{analind}.parentdir = NewRootDir;
            NM.analysis{analind}.rootdir = NewAnalDir;
            NM.analysis{analind}.logfile = fullfile(NewAnalDir,['NM_Analysis' NM.analysis{analind}.id '.log']);
            [~,paramsold] = fileparts(NM.analysis{analind}.paramfile);
            NM.analysis{analind}.paramdir = fullfile(NewAnalDir,'params');
            NM.analysis{analind}.paramfile = fullfile(NewAnalDir,'params',[paramsold '.mat']);

            % Adjust CVdatamat directories
            algostr = NM.analysis{analind}.params.TrainParam.SVM.prog;
            for i=1:numel(NM.analysis{analind}.GDdims)
                NM.analysis{analind}.GDdims{i}.RootPath = fullfile(NewAnalDir,algostr);
            end

            % Adjust OOCVdatamat directories
            if isfield(NM.analysis{analind},'OOCV')
                for i=1:numel(NM.analysis{analind}.OOCV)
                    oocvdir = sprintf('OOCV_%g',i);
                    NM.analysis{analind}.OOCV{i}.RootPath = fullfile(NewAnalDir,algostr,oocvdir);
                end
            end
         end
         succeeded{cnts} = NM.analysis{analind}.id;
         cnts=cnts+1;
     catch
         failed{cntf} = NM.analysis{analind}.id;
         cntf=cntf+1;
     end
end
succeededstr = ''; failedstr='';
if ~isempty(succeeded)
    succeededstr = sprintf('%s',strjoin(succeeded,'\n'));
    succeededstr = sprintf('The following analyses have been successfully adjusted to the new root directory:\n%s', succeededstr);
end
if ~isempty(failed)
    failedstr = sprintf('%s',strjoin(failed,'\n'));
    failedstr = sprintf('\n\nThe following analyses could not be adjusted to the new root directory:%s', failedstr);
end
msgbox(sprintf('%s%s',succeededstr,failedstr));

