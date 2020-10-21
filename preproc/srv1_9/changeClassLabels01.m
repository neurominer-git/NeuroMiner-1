function formatedClassLabels=changeClassLabels01(originalClassLabels)
% Change class labels to 0,1,2,...
% Usage:
% formatedClassLabels=changeClassLabels01(originalClassLabels)
% originalClassLabels: column vector of numbers of string cells.
% formatedClassLabels: column vector including the 0,1,2,3,...class labels 
% for example,
% [-1;-1;1;1] to [0;0;1;1]; [-1;-1;1;1;2;2] to [0;0;1;1;2;2]
% {'normal';'normal';'cancer';'cancer'} to [0;0;1;1]
% Contact Information:
% Yifeng Li
% University of Windsor
% li11112c@uwindsor.ca; yifeng.li.cn@gmail.com
% May 18, 2011

formatedClassLabels=nan(size(originalClassLabels));
uniqueLabels=unique(originalClassLabels);
for c=1:numel(uniqueLabels)
    if isnumeric(uniqueLabels)
        formatedClassLabels(originalClassLabels==uniqueLabels(c))=c-1;
    end
    if iscellstr(uniqueLabels)
        formatedClassLabels(strcmp(originalClassLabels,uniqueLabels{c}))=c-1;
    end
end

end