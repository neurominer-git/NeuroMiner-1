function [X, L] = nk_CompatY2(NM, varind, oocvind, FUSION)
% =========================================================================
% FORMAT [X, L] = nk_CompatY2(NM, varind, oocvind)
% =========================================================================
% This function extracts the feature space from the NM structure, which is 
% needed for preprocessing, training and visualization. The extraction is 
% controlled by the FUSION settings defined by the user. 
% If OOCV data is present and an OOCV index is given the returned 
% feature space container will contain also the respectively processed OOCV 
% dataset.
%
% Inputs:
% -------
% NM                : the NM structure
% varind            : Index vector to modalities
% oocvind           : Index to OOCV data container
% FUSION            : Data fusion parameters
%
% Outputs:
% --------
% X                 : Either a single struct or a struct array depending on
%                     whether FUSION.flag = 0|1|3 or 2 (intermediate fusion)
% X.dimvecx         : Vector of starting column of each concatenated
%                     modality
% X.dimsizes        : Feature space sizes of each modality
% X.fulldimsizes    : Full feature space sizes of each modality (before
%                     feature extraction using badcoords vector)
% X.Y               : Feature space(s) depending on FUSION.M and
%                     FUSION.flag
% X.brainmask       : Brainmask for each modality (if modality is
%                     neuroimaging data)
% X.badcoords       : Vector of zero/low var voxels for each modality 
%                     (if modality is neuroimaging data)
% X.datadesc        : Modality descriptor
% X.datatype        : Data type descriptor
% X.threshval       : Threshold (if modality is neuroimaging data)
% X.threshop        : Threshold operation (if modality is neuroimaging data)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris 10/2020

if ~exist('varind','var'), varind=1; end; nM = numel(varind);
if ~exist('oocvind','var')
   if isfield(NM,'OOCV')
       oocvflag = nk_input('Do you want to extract OOCV data?',0,'yes|no',[1,0],2);
       if oocvflag
           [~,~,oocvind] = nk_SelectOOCVdata(NM,1,0);
           if isempty(oocvind), warning('No OOCV data selected. Skipping OOCV extraction!'); end 
       end
   else
       oocvind=[]; 
   end
end
if ~exist('FUSION','var') || isempty(FUSION) 
    if nM > 1
        [~, FUSION.flag] = nk_Fusion_config(NM, varind);
    else
        FUSION.flag = 0;
    end
    FUSION.M = varind;
end
L = NM.label;

%Load linked OOCV files if needed
if ~isempty(oocvind) && ischar(NM.OOCV{oocvind}.Y) 
   if ~exist(NM.OOCV{oocvind}.Y,'file')
       error('Linked OOCV container not found in expected path: %s',NM.OOCV{oocvind}.Y);
   end
   fprintf('\nLoading linked OOCV container: %s',NM.OOCV{oocvind}.Y);
   load(NM.OOCV{oocvind}.Y);
   NM.OOCV{oocvind} = OOCV;
end

switch FUSION.flag

    case {0,3}
        
        X.dimvecx(1) = 0;
        X.Y = NM.Y{varind};
        X.brainmask{1} = NM.brainmask{varind};
        if ~isempty(X.brainmask{1})
            X.brainmask{1} = nk_RemCommaFromStr( X.brainmask{1} );
        end
        if ~isempty(oocvind) && isfield(NM,'OOCV') && oocvind <= numel(NM.OOCV),
            X.Yocv = NM.OOCV{oocvind}.Y{varind}; 
            % If the original space-defining image cannot be found use
            % the one used to read-in the OOCV data
            if ~isempty(X.brainmask{1}) && ~strcmp(X.brainmask{1},'none') && ~exist(X.brainmask{1},'file')
                X.brainmask{1} = nk_RemCommaFromStr( NM.OOCV{oocvind}.datadescriptor{varind}.input_settings.brainmask );
            end
        end
        if ~isempty(oocvind) && isfield(NM,'C')  && oocvind <= numel(NM.C)
            X.C = NM.C{oocvind}.Y{varind}  ;
            if isfield(NM.C{oocvind},'label'), X.Clabels = NM.C{oocvind}.label; end
        end
        X.dimsizes  = size(X.Y,2);
        X.datadesc = NM.datadescriptor{varind}.desc;
        X.datatype = NM.datadescriptor{varind}.type;
        if isfield(NM.datadescriptor{varind},'Vm')
            X.Vm = NM.datadescriptor{varind}.Vm;
            X.Vmvol = NM.datadescriptor{varind}.Vmvol;
        end
        if isfield(NM.datadescriptor{varind},'threshval')
            X.threshval = NM.datadescriptor{varind}.threshval;
            X.threshop = NM.datadescriptor{varind}.threshop{1};
            if isfield( NM.datadescriptor{varind},'desc_label' ), X.threshlb = NM.datadescriptor{varind}.desc_label; end
        else
            X.threshval = 0;X.threshop = 'gt';
        end
        if isfield( NM.datadescriptor{varind},'Yw' ), X.Yw = NM.datadescriptor{varind}.Yw; end
        X.badcoords{1} = NM.badcoords{varind};
        if isempty(X.badcoords{1}), X.badcoords{1} = zeros(1,X.dimsizes); end
        X.dimvecx(2) = X.dimsizes; 
         
    case 1 % Concatenate modalities
        
        X.dimvecx   = 0;
        X.Y = [];
        
        for i = 1:nM

            % Training and CV data
            X.Y = [ X.Y NM.Y{varind(i)} ] ;
            X.brainmask{i} = nk_RemCommaFromStr( NM.brainmask{varind(i)} );
            % Independent test data
            if ~isempty(oocvind) && isfield(NM,'OOCV')  && oocvind <= numel(NM.OOCV)
                X.Yocv = [ X.Yocv NM.OOCV{oocvind}.Y{varind(i)} ] ;
                % If the original space-defining image cannot be found use
                % the one used to read-in the OOCV data
                if ~isempty(X.brainmask{i}) && ~strcmp(X.brainmask{i},'none') && ~exist(X.brainmask{i},'file')
                    X.brainmask{i} = nk_RemCommaFromStr( NM.OOCV{oocvind}.datadescriptor{varind(i)}.input_settings.brainmask );
                end
            end
            % Calibration data
            if ~isempty(oocvind) && isfield(NM,'C')  && oocvind <= numel(NM.C)
                X.C = [ X.C NM.C{oocvind}.Y{varind(i)} ] ;
                if isfield(NM.C{oocvind},'label'), X.Clabels = NM.C{oocvind}.label; end
            end
            X.datatype(i) = NM.datadescriptor{varind(i)}.type;
            X.datadesc{i} = NM.datadescriptor{varind(i)}.desc;
            if isfield(NM.datadescriptor{varind(i)},'Vm')
                X.Vm{i} = NM.datadescriptor{varind(i)}.Vm;
                X.Vmvol{i} = NM.datadescriptor{varind(i)}.Vmvol;
            end
            if isfield(NM.datadescriptor{varind(i)},'threshval')
                X.threshval{i} = NM.datadescriptor{varind(i)}.threshval;
                X.threshop{i} = NM.datadescriptor{varind(i)}.threshop{1};
                if isfield( NM.datadescriptor{varind(i)},'desc_label' ) && ...
                        ~isempty( NM.datadescriptor{varind(i)}.desc_label )
                    X.threshlb{i} = NM.datadescriptor{varind(i)}.desc_label; 
                end
            else
                X.threshval{i} = 0; X.threshop{i} = 'gt';
            end
            X.dimsizes(i) = size(NM.Y{varind(i)},2);
            if isfield( NM.datadescriptor{varind(i)},'Yw' ), X.Yw{i} = NM.datadescriptor{varind(i)}.Yw; end
            X.badcoords{i} = NM.badcoords{varind(i)};
            if isempty(X.badcoords{i}), X.badcoords{i} = zeros(1,X.dimsizes(i)); end
            X.fulldimsizes(i) = size(NM.badcoords{varind(i)},2);
            X.dimvecx(i+1) = sum(X.dimsizes(1:i)); 
        
        end

    case 2 % store modalities in cell array
        
        for i=1:nM

            X(i).dimsizes  = size(NM.Y{varind(i)},2);
            X(i).fulldimsizes(i) = size(NM.badcoords{varind(i)},2);
            X(i).dimvecx   = [ 0 X(i).dimsizes ];
            X(i).Y         = NM.Y{varind(i)};
            X(i).brainmask{1} = nk_RemCommaFromStr( NM.brainmask{varind(i)} );
            if ~isempty(oocvind) && isfield(NM,'OOCV')  && oocvind <= numel(NM.OOCV)
                X(i).Yocv = NM.OOCV{oocvind}.Y{varind(i)};
                % If the original space-defining image cannot be found use
                % the one used to read-in the OOCV data
                if ~isempty(X(i).brainmask{1}) && ~strcmp(X(i).brainmask{1},'none') && ~exist(X(i).brainmask{1},'file')
                    X(i).brainmask{1} = nk_RemCommaFromStr( NM.OOCV{oocvind}.datadescriptor{varind(i)}.input_settings.brainmask );
                end
            end
            if ~isempty(oocvind) && isfield(NM,'C')  && oocvind <= numel(NM.C)
                X(i).C = NM.C{oocvind}.Y{varind(i)} ;
                if isfield(NM.C{oocvind},'label'), X(i).Clabels = NM.C{oocvind}.label; end
            end
           
            X(i).datatype = NM.datadescriptor{varind(i)}.type;
            X(i).datadesc = NM.datadescriptor{varind(i)}.desc;
            if isfield(NM.datadescriptor{varind(i)},'Vm')
                X(i).Vm = NM.datadescriptor{varind(i)}.Vm;
                X(i).Vmvol = NM.datadescriptor{varind(i)}.Vmvol;
            end
            if isfield(NM.datadescriptor{varind(i)},'threshval')
                X(i).threshval = NM.datadescriptor{varind(i)}.threshval;
                X(i).threshop = NM.datadescriptor{varind(i)}.threshop{1};
                if isfield( NM.datadescriptor{varind(i)},'desc_label' ) && ...
                        ~isempty( NM.datadescriptor{varind(i)}.desc_label ), 
                    X(i).threshlb = NM.datadescriptor{varind(i)}.desc_label; 
                end
            else
                X(i).threshval = 0; X(i).threshop = 'gt';
            end
            if isfield( NM.datadescriptor{varind(i)},'Yw' ), X(i).Yw = NM.datadescriptor{varind(i)}.Yw; end
            X(i).badcoords{1} = NM.badcoords{varind(i)};
            if isempty(X(i).badcoords{1}), X(i).badcoords{1} = zeros(1,X(i).dimsizes); end
        end
end