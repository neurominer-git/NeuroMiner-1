% script to test the adjustment to combat

%matlab
p=10000;
n=10;
batch = [1 1 1 1 1 2 2 2 2 2]; %Batch variable for the scanner id
dat = randn(p,n); %Random data matrix

%and let simulate an age and disease variable: 
age = [82 70 68 66 80 69 72 76 74 80]'; % Continuous variable
disease = [1 2 1 2 1 2 1 2 1 2]'; % Categorical variable

%We create a n x 2 model matrix with age as the first column, and the second disease group as a dummy variable for the second column (disease=1 being the baseline category):
disease = dummyvar(disease);
mod = [age disease(:,2)];

%We use the function `combat` to harmonize the data across the 2 scanners:
[data_harmonized,estimators] = combat(dat, batch, mod);

%Test the same data with the estimators 
[data_harmonized2,estimators] = combat(dat, batch, mod, estimators); 
isequal(data_harmonized,data_harmonized2) % equal

% Test with a reduced sample and compare to the original sample correction.
dat2 = dat(:,1:5);
batch2 = batch(1:5);
mod2 = mod(1:5,:); 
[data_harmonized3,estimators] = combat(dat2, batch2, mod2, estimators); 
isequal(data_harmonized(:,1:5),data_harmonized3(:,1:5)) % it is equal
