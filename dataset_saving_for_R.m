%%%%% SAVE NEW DATASET WITH NANs for plotting in R

%% Load data
data = readtable("Patient_Master.csv");

%% Dataset sub-divisions
idx_samples.HC = find(string(data.COHORT)=='HC'); % Healthy Controls
idx_samples.PD = find(string(data.COHORT)=='PD'); % Parkinson's Disease
idx_samples.SWEDD = find(string(data.COHORT)=='SWEDD'); % Scans without evidence of dopaminergic deficit
idx_samples.Prodromal = find(string(data.COHORT)=='Prodromal'); % Early signs of PD - NO ONE

data_pd = data([idx_samples.Prodromal; idx_samples.PD; idx_samples.SWEDD],:);
data_hc = data(idx_samples.HC, :);

%% Saving new datasets with only principal variables:
% np1ptot, np1rtot, np2tot, np3tot, np4tot, genetics, familiarity, ethnicity, sex, age,
% height, weight, hand, primary diagnosis
new_data_hc = table();
new_data_pd = table();

important_var = [34, 41, 55, 94, 104, 5, 11, 7, 8, 4, 153, 154, 9, 12];
names = data_hc.Properties.VariableNames;

for i = 1:length(important_var)
    new_data_hc(:,i) = data_hc(:,names(important_var(i)));
    new_data_pd(:,i) = data_pd(:,names(important_var(i)));
end

new_data_hc.Properties.VariableNames = names(1, important_var);
new_data_pd.Properties.VariableNames = names(1, important_var);

writetable(new_data_hc, 'new_data_hc.csv');
writetable(new_data_pd, 'new_data_pd.csv');

clear new_data_pd new_data_hc