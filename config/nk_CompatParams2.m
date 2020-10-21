function [PREPROC, ...
            RFE, ...
            GRD, ...
            SVM, ...
            LIBSVMTRAIN, ...
            LIBSVMPREDICT, ...
            RVM, ...
            MKLRVM, ...
            CMDSTR, ...
            MULTI, ...
            VIS, paramstr] = nk_CompatParams2(TrainParam, varind, paramstr)

        global MODEFL
        
PREPROC             = [];
RFE                 = []; 
GRD                 = [];
SVM                 = [];
LIBSVMTRAIN         = [];
LIBSVMPREDICT       = [];
RVM                 = [];
MKLRVM              = [];
CMDSTR              = [];
VIS                 = [];
MULTI               = [];
if ~exist('paramstr','var'), paramstr = []; end

if ~isempty(TrainParam.FUSION) && TrainParam.FUSION.flag == 3
   TrainParam = TrainParam.STRAT{varind};
end

if isfield(TrainParam,'PREPROC'),
    if iscell(TrainParam.PREPROC) 
        if varind > numel(TrainParam.PREPROC);
            warning('VARIND out of bounds. Resetting to VARIND = 1 !!!')
            varind = 1;
        end
        if numel(varind)>1
            PREPROC  = TrainParam.PREPROC(varind);
        else
            PREPROC  = TrainParam.PREPROC{varind}; 
        end
    else
        PREPROC  = TrainParam.PREPROC; 
    end
end

if isfield(TrainParam,'VIS'),
    if iscell(TrainParam.VIS) 
        if varind > numel(TrainParam.VIS);
            warning('VARIND out of bounds. Resetting to VARIND = 1 !!!')
            varind = 1;
        end
        if numel(varind)>1
            VIS= TrainParam.VIS(varind);
        else
            VIS = TrainParam.VIS{varind}; 
        end
    else
        VIS = TrainParam.VIS; 
    end
end


if isfield(TrainParam,'SVM'), SVM = TrainParam.SVM; end
SVM.RVMflag = 0;
            
switch SVM.prog
    case 'LIBSVM'
        try
            CMDSTR  = nk_DefineCmdStr(SVM, MODEFL);
            if SVM.LIBSVM.Optimization.b, SVM.RVMflag = true; end
            switch SVM.LIBSVM.LIBSVMver
                case 0
                    LIBSVMTRAIN = '312'; LIBSVMPREDICT = '312';
                case 2
                    LIBSVMTRAIN = '289'; LIBSVMPREDICT = '289';
                case 1
                    LIBSVMTRAIN = '291'; LIBSVMPREDICT = '291';
                case 3
                    LIBSVMTRAIN = '289PLUS'; LIBSVMPREDICT = '289PLUS';
            end
            LIBSVMTRAIN = ['svmtrain' LIBSVMTRAIN]; LIBSVMPREDICT = ['svmpredict' LIBSVMPREDICT];
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Parameters for LIBSVM');
        end
        
    case 'LIBLIN'
        try
            CMDSTR  = nk_DefineCmdStr(SVM, MODEFL);
            if SVM.LIBLIN.b, SVM.RVMflag = true; end
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Parameters for LIBLINEAR');
        end
        
    case 'MikRVM'
        try
            RVM.UserOpt = TrainParam.SVM.RVM.UserOpt;
            RVM.ParamSet = TrainParam.SVM.RVM.ParamSet;
            RVM.LikelihoodModel = TrainParam.SVM.RVM.LikelihoodModel;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Parameters for Mike Tipping RVM');
        end
        SVM.RVMflag = true;
    case 'MKLRVM'
        try
            MKLRVM  = TrainParam.SVM.MKLRVM;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Multiple Kernel Learning parameters for RVM');
        end
        SVM.RVMflag = 1;
    case {'BLOREG', 'IMRELF', 'kNNMEX','GLMFIT','GLMNET'}
        SVM.RVMflag = 1;
    case 'RNDFOR'
        CMDSTR  = nk_DefineCmdStr(SVM, MODEFL);
        SVM.RVMflag = 1;
    otherwise
        CMDSTR  = nk_DefineCmdStr(SVM);

end

if isfield(TrainParam,'GRD'),            GRD      = TrainParam.GRD; end
if isfield(TrainParam,'RFE'),            RFE      = TrainParam.RFE; end
if isfield(TrainParam,'MULTI'),          MULTI    = TrainParam.MULTI; end