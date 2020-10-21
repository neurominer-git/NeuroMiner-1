function [] = calculateComponentWeightedAverageNIFTI(dataList,resultsDir,numBases)

% dataList: .csv file containing the images for which the coefficients need
%           to be calculated. Full path for every image file is given in
%           every line
% resultsDir: directory where the NMF results are (i.e., the level where 
%             the NumBases folder is placed)
% numBases: determines the solution for which one wants to calculate 
%           subject coefficients

for b=1:length(numBases)
    disp(numBases(b))
    dataPath=[resultsDir '/NumBases' num2str(numBases(b)) '/OPNMF/niiImg/'];        
    
    % loading estimated non-negative components
    listing = dir(dataPath);
    hh =cellfun(@(x) ( (strfind(x,'Basis')==1)  ),{listing(:).name},'UniformOutput',false) ;
    listing=listing(cellfun(@(x) (~isempty(x)&&(x==1)),hh));
    
    if(length(listing)~=numBases(b))
        error(['I cannot find ' num2str(numBases(b)) ' basis images in the folder ' dataPath ])
    end
    
    for i=1:numBases(b)
        nii = load_untouch_nii([dataPath listing(i).name]);
        B(:,i) = double(nii.img(:)');        
    end
    
    % normalize to sum to 1
    Blen = sum(B,1);
    if any(Blen==0)
        Blen(Blen==0) = 1;
    end
    nB = bsxfun(@times,B,1./Blen) ;
        
    % since the size and number of files is such that we can not manage
    % in batch mode, we are going to calculate weighted average values
    % subject by subject
    
    % read list
    fid=fopen(dataList,'r');
    if (fid == -1)
        error(['extractBases:calculateComponentWeightedAverage ','Can not open ' list ' file.']);
    end
    datafullpath = textscan(fid,'%s\n');
    fclose(fid);
    
    datafullpath = datafullpath{1,1} ;
    datafullpath = cellstr(datafullpath) ;
    count = numel(datafullpath);
    
    fid = fopen([resultsDir  '/NumBases' num2str(numBases(b)) '/cmpWeightedAverageNumBases_' num2str(numBases(b)) '.csv'],'w');
    frmtWrite='%s,';
    frmtWrite=[frmtWrite repmat('%f,',1,numBases(b)-1)];
    frmtWrite=[frmtWrite '%f\n'];
    
    wA = zeros(count,numBases) ;
    for i=1:count
        nii = load_untouch_nii(datafullpath{i});
        wA(i,:) = double(nii.img(:)')*nB;
        fprintf(fid,frmtWrite,datafullpath{i},wA(i,:)');
    end
    fclose(fid);         
end