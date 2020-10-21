function [P, Y] = nk_FSreslice(S,P)

Y= [];

if ~exist('S','var') || isempty(S)
    S = spm_select(Inf,'dir','Specify directories/list of subjects to work on');
    [S,SDIR] = spm_str_manip(S,'c');
    if ~exist(fullfile(SDIR,deblank(S(1,:))),'dir')
        SDIR = spm_select(1,'dir','Specify FREESURFER subject directory');
    end
end

if ~exist('P','var') || isempty(P)
    
    P.S = S;
    
    P.SDIR = SDIR;
    
    P.FSPATH = spm_select(1,'dir','Freesurfer main directory');
    
    P.HEMI = nk_input('Hemisphere',0,'rh|lh|both',[1,2,3],1);

    P.TARGET = nk_input('Target space',0,'m', ...
                        ['fsaverage4|' ...
                          'fsaverage5|' ...
                          'fsaverage6|'],1:3,1);

    P.SURFVAL = nk_input('Value of measure',0,'m', ...
                        ['thickness|', ...
                         'area|', ...
                         'area.pial|', ...
                         'sulc|', ...
                         'curv|', ...
                         'jacobian.white|', ...
                         'volume'],1:7,1);

    P.SURFTYPE = nk_input('Type of measure',0,'m', ...
                          ['Freesurfer curvature file|', ...
                           'Freesurfer paint file|', ...
                           'other (as accpeted by mri_convert'],1:3,1);


    P.FWHM_IN = nk_input('Source space Gaussian smoothing',0,'e');

    P.FWHM_OUT = nk_input('Target space Gaussian smoothing',0,'e');

    P.OUTTYPE = nk_input('Output data type',0,'m', ...
                         ['mgz|' ...
                          'mgh|' ...
                          'nifti|' ...
                          'analyze'], 1:4); 

    P.TRGIORDER = nk_input('Icosahedron order',0,'e',5);

    P.OUTPATH = spm_select(1,'dir','Output directory');

end

switch P.HEMI
    case 3
        hemi = {'rh','lh'};
    case 1
        hemi = {'rh'};
    otherwise 
        hemi = {'lh'};
end

trgsubject = {'fsaverage4','fsaverage5','fsaverage6'};
trgsubject = trgsubject{P.TARGET};

surfval = {'thickness','area','area.pial','sulc','curv','jacobian.white','volume'};
surfval = surfval{P.SURFVAL};

surftype = {'curv','w'};
surftype = surftype{P.SURFTYPE};

fwhm_in = sprintf('%g',P.FWHM_IN);
fwhm_out = sprintf('%g',P.FWHM_OUT);

outtype = {'mgz','mgh','nifti3d','analyze3d'};
outext = {'mgz','mgh','img','img'};
outtype = outtype{P.OUTTYPE};
outext = outext{P.OUTTYPE};
trgiorder = sprintf('%g',P.TRGIORDER);

%% CREATE FREESURFER SCRIPT FOR SURFACE2SURFACE REGISTRATION AND DOWNSAMPLING
fid = fopen('Surf2Surf.sh','w');

fprintf(fid,'export SUBJECTS_DIR=%s',P.SDIR);
fprintf(fid,'\nexport FREESURFER_HOME=%s',P.FSPATH);
fprintf(fid,'\nsource $FREESURFER_HOME/SetUpFreeSurfer.sh');
cnt=1;
for h = 1:numel(hemi)

    for i=1:size(S,1)
        Si = regexprep(deblank(S(i,:)),'/','');
        trgfile{cnt} = sprintf('%s-%s-%s-%s-fwhmin%s-fwhmout%s.%s',Si,hemi{h},surfval,trgsubject,fwhm_in,fwhm_out,outext);
        fprintf(fid,['\n mri_surf2surf', ...
                        ' --hemi %s', ...
                        ' --srcsubject %s', ...
                        ' --srcsurfval %s', ...
                        ' --src_type %s', ...
                        ' --trgsubject %s', ...
                        ' --trgicoorder %s', ...
                        ' --trgsurfval %s', ...
                        ' --trg_type %s', ...
                        ' --fwhm-src %s', ...
                        ' --fwhm-trg %s'], ...
                        hemi{h}, ...
                        Si, ...
                        surfval, ...
                        surftype, ...
                        trgsubject, ...
                        trgiorder, ...
                        trgfile{cnt}, ...
                        outtype, ...
                        fwhm_in, ...
                        fwhm_out);
        cnt = cnt+1;

    end

end

%% SCRIPT EXECUTION
execflag = nk_input('Execute script?',0,'yes|no',[1,0],1);

if execflag
    system('sh Surf2Surf.sh');
end

file2matflag = nk_input('Read in files and construct data matrix?',0,'yes|no',[1,0],1);
if file2matflag
    cnt = 1;
    H = cell(numel(hemi),1);
    for h = 1:numel(hemi)
        for i=1:size(S,1)
            filename = fullfile(pwd,trgfile{cnt});
            fprintf('\nReading file: %s',filename)
            switch outtype
                case {'mgz','mgh'}
                    y = MRIread(filename);
                    H{h}=[H{h}; y.vol];
                case {'nifti3d','nifti4d'}
                    y = spm_read_vols(spm_vol(filename));
                    H{h}=[H{h}; y'];
            end
            cnt=cnt+1;
        end
    end
    for h = 1:numel(hemi)
        Y = [Y H{h}];
    end
end
end

    