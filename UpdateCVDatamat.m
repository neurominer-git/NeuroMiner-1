function UpdateCVDatamat

spm5ver = nk_CheckSPMver;

if spm5ver 
    P = spm_select(Inf,'.*CVdatamat.*mat','Select CVdatamats');
else
    P = spm_get(Inf,'*CVdatamat.*mat','Select CVdatamats');
end

for i=1:size(P,1)
   pth = deblank(P(i,:));
   fprintf('\nLoading: %s',pth)
   load(pth)
   if ~exist('GD','var')
       fprintf('\nThis is not a valid CVdatamat.')
   else
        for h=1:size(GD,1)
            if isfield(GD{h},'RT'), GD{h} = rmfield(GD{h},'RT'); end
            if isfield(GD{h},'RV'), GD{h} = rmfield(GD{h},'RV'); end
            if isfield(GD{h},'RS'), GD{h} = rmfield(GD{h},'RS'); end
            if isfield(GD{h},'Mm'), GD{h} = rmfield(GD{h},'Mm'); end
            if isfield(GD{h},'Md'), GD{h} = rmfield(GD{h},'Md'); end
            if isfield(GD{h},'W2'), GD{h} = rmfield(GD{h},'W2'); end
            if isfield(GD{h},'mMd'), GD{h} = rmfield(GD{h},'mMd'); end
            if isfield(GD{h},'bestw2'), GD{h} = rmfield(GD{h},'bestw2'); end
            if isfield(GD{h},'bestmarg'), GD{h} = rmfield(GD{h},'bestmarg'); end
            if isfield(GD{h},'bestnorm'), GD{h} = rmfield(GD{h},'bestnorm'); end
%             GD{h}.Mm = []; 
%             GD{h}.Md = []; 
%             GD{h}.W2=[]; 
%             GD{h}.mMd = [];
%             GD{h}.bestw2 = []; 
%             GD{h}.bestmarg = [];
%             GD{h}.bestnormw = [];
        end
        fprintf('\nSaving: %s',pth)
        if exist('MD','var')
            save(pth,'GD','MD','k','g','operm','ofold');
        else
            save(pth,'GD','k','g','operm','ofold');
        end
   end
   clear GD MD k g operm ofold
end
end