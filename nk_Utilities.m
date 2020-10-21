function nk_Utilities
global CV NM

menustr = ['Set root paths of neuroimaging tools (SPM/Freesurfer)|' ...
           'Create PreprocData path master|' ...
           'Create CVdatamat path master|' ...
           'Create CVresults path master|' ...
           'Create VISdatamat path master|' ...
           'Create OptPreprocParam path master|' ...
           'Create OptModel path master|' ...
           'Update analyses'' paths to new root directory'];
menuact = 2:9;

nk_PrintLogo
act = nk_input('Choose utility function',0,'mq', menustr, menuact, 1);
                 
switch act
    case 0
        return
    case 2
        neurominerpath = fileparts(which('neurominer.m'));
        imaging_init_path = fullfile(neurominerpath,'imaging_init.mat');
        nk_ImagingInit(neurominerpath, imaging_init_path, 1);
    case 3
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenPreprocMaster2(NM.id, CV(1), [], inp.rootdir, [], 1, inp.varind, inp.varstr, inp.concatfl);
    case 4
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenCVdataMaster2(NM.id, CV(1), [], inp.rootdir, [], 1, inp.varind, inp.varstr, inp.concatfl);
    case 5
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenCVresultsMaster(NM.id,[],inp.rootdir,[],1); 
    case 6
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenDataMaster(NM.id, 'VISdatamat', CV(1),[],inp.rootdir,[],1);
    case 7
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenDataMaster(NM.id, 'OptPreprocParam', CV(1),[],inp.rootdir,[],1);
    case 8
        inp = []; inp = nk_GetAnalModalInfo_config(NM, inp);
        nk_GenDataMaster(NM.id, 'OptModel', CV(1),[],inp.rootdir,[],1);
    case 9
        complvec = []; for z=1:numel(NM.analysis), if NM.analysis{z}.status, complvec = [ complvec z ]; end; end
        t_act = 1; brief = 1; analind = 1; showmodalvec = []; 
        while t_act>0, 
            [t_act, analind, ~, showmodalvec , brief] = nk_SelectAnalysis(NM, 0, 'MAIN INTERFACE >> UPDATE ANALYSES ROOT DIRECTORIES ', analind, [], 1, showmodalvec, brief); 
        end
        if ~isempty(analind), 
            analind = complvec(analind);
        else
            analind = complvec;
        end
        newdir = nk_DirSelector('Update analyses'' root paths');
        NM = nk_UpdateRootPaths(NM, analind, newdir);
    case 10
        [P, Y] = nk_FSreslice;
        assignin('base','P',P);
        assignin('base','Y',Y);    
end

nk_Utilities

end

