% function that evaluates the reconstruction error only on whole brain
clear
clc 

% read the data
data_lst = '/cbica/home/chenyin/PD_Spatial_Project/Script/Atlas_CER_W/CER_W_All.csv';
fid=fopen(data_lst,'r');
datafullpath = textscan(fid,'%s %d\n');
fclose(fid);

datafullpath = datafullpath{1,1} ;

count = numel(datafullpath);
info = load_untouch_header_only(datafullpath{1});

% load data
sizeImg = info.dime.dim(2)*info.dime.dim(3)*info.dime.dim(4) ;
data.X = zeros(sizeImg,count);
for i=1:count
    disp(i/count)
    nii = load_untouch_nii(datafullpath{i});    
    data.X(:,i) = nii.img(:) ;
end

% load results and calculate reconstruction error
% path were the resutls have been stored
resultsPath='/cbica/home/chenyin/PD_Spatial_Project/Script/Atlas_CER_W/sge_job_output/';

% number of components for which solutions have been calculated
numBases=2:2:30;

RecError=zeros(length(numBases),1);
for b=1:length(numBases)
    disp(b/length(numBases))
    load([resultsPath 'NumBases' num2str(numBases(b)) '/OPNMF/ResultsExtractBases.mat'])  
    Est = B*C ;
    RecError(b) = norm(data.X-Est,'fro') ;    
    clear B C
end

% make figure
outputDir='/volume/GI_HCHR_RS/Data/NMF_FS/test/split1/';

% 1) reconstruction error
figure;plot(numBases,RecError,'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Reconstruction error','fontsize',12)
xlim([numBases(1) numBases(end)])
set(gca,'fontsize',12)
saveas(gcf,[outputDir 'RecError.fig'])
saveas(gcf,[outputDir 'RecError.png'])

% 2) gradient of reconstruction error
figure;plot(numBases(2:end),diff(RecError),'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Gradient of reconstruction error','fontsize',12)
xlim([numBases(1) numBases(end)])
set(gca,'fontsize',12)
saveas(gcf,[outputDir 'gradientRecError.fig'])
saveas(gcf,[outputDir 'gradientRecError.png'])

% 3) Percentage of improvement over range of components used
figure;plot(numBases,abs(RecError-RecError(1))./abs(RecError(1)-RecError(end)),'r','LineWidth',2)
xlabel('Number of components','fontsize',12)
ylabel('Percentage of improvement over range of components used','fontsize',12)
xlim([numBases(1) numBases(end)])
set(gca,'fontsize',12)
saveas(gcf,[outputDir 'percentageImprovementRecError.fig'])
saveas(gcf,[outputDir 'percentageImprovementRecError.png'])

close all
