clear all
close all
clc

%% load data
data = readtable("Patient_Master.csv");

DAT_SCAN_SPECT = [data.DATSCAN_CAUDATE_R data.DATSCAN_CAUDATE_L data.DATSCAN_PUTAMEN_R  data.DATSCAN_PUTAMEN_L data.DATSCAN_PUTAMEN_R_ANT data.DATSCAN_PUTAMEN_L_ANT];
%DAT_SCAN_PET =  [data.AV133_RCAUD_S data.AV133_LCAUD_S data.AV133_RPUTANT_S data.AV133_RPUTPOST_S data.AV133_LPUTANT_S data.AV133_LPUTPOST_S];

ROIs_labels = ["Right Caudate", "Left Caudate", "Right Putamen", "Left Putamen", "Left Anterior Putamen", "Right Anterior Putamen"];

idx_HC = find(string(data.COHORT)=='HC');
HC_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_HC,:);

idx_PD = find(string(data.COHORT)=='PD');
PD_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_PD,:);

idx_SWEDD = find(string(data.COHORT)=='SWEDD');
SWEDD_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_SWEDD,:);

idx_Prodromal = find(string(data.COHORT)=='Prodromal');
Prodromal_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_Prodromal,:);

NON_HC_DAT_SCAN_SPECT = DAT_SCAN_SPECT([idx_Prodromal; idx_PD; idx_SWEDD],:);

%% histograms
for i=1:6
    subplot(2,3,i)
    histogram(NON_HC_DAT_SCAN_SPECT(:,i),'FaceColor','auto','Normalization','probability')
    xlabel('Striatal binding ratio [adim]')
    title(['SBR in ' ROIs_labels(i)])
    hold on 
    mu_non_hc = mean(NON_HC_DAT_SCAN_SPECT(:,i));
    xline(mu_non_hc,'LineWidth',2,'Color','r')
    histogram(HC_DAT_SCAN_SPECT(:,i),'FaceColor','auto','Normalization','probability')
    mu_hc = mean(HC_DAT_SCAN_SPECT(:,i));
    xline(mu_hc,'LineWidth',2,'Color','b')
    xlabel('Striatal binding ratio [adim]')
    title(['SBR in ' ROIs_labels(i)])
    hold off
    legend('Non healthy','Average DAT non healthy','Healthy controls','Average DAT healthy')
end

%% Analysis missing data ex family history
answers_fam_pd = data.ANYFAMPD;
idx_yes_fam_pd = find(string(answers_fam_pd) =='1');
%yes_fam_pd = answers_fam_pd
num_yes_fam_pd = length(idx_yes_fam_pd);
idx_no_fam_pd = find(string(answers_fam_pd) =='0');
num_no_fam_pd = length(idx_no_fam_pd);
idx_nan_fam_pd = find(string(answers_fam_pd) =='NA');
num_nan_fam_pd = length(idx_nan_fam_pd);

percent_not_nan = (num_no_fam_pd + num_yes_fam_pd)/(num_no_fam_pd + num_yes_fam_pd+num_nan_fam_pd);

%% Vari riassuntive
np1r = data.NP1RTOT;
np1p = data.NP1PTOT;
np2 = data.NP2PTOT;
np3 = data.NP3TOT;
np4 = data.NP4TOT;

% Find missing data
idx_nan_np1r = find(isnan(np1r));
idx_nan_np1p = find(isnan(np1p));
idx_nan_np2 = find(isnan(np2));
idx_nan_np3 = find(isnan(np3));
idx_nan_np4 = find(string(np4) =='NA'); 
idx_nan_np1 = union(idx_nan_np1r, idx_nan_np1p);
idx_nan_np23 = union(idx_nan_np2, idx_nan_np3);
idx_nan_np123 = union(idx_nan_np23,idx_nan_np1);
idx_nan_np = union(idx_nan_np123,idx_nan_np4);

data_without_nan_np = data(:,:);
data_without_nan_np(idx_nan_np123,:) = [];

np1r = data_without_nan_np.NP1RTOT;
np1p = data_without_nan_np.NP1PTOT;
np2 = data_without_nan_np.NP2PTOT;
np3 = data_without_nan_np.NP3TOT;
np4 = data_without_nan_np.NP4TOT;

% Normalization
np1r_norm = np1r./sum(np1r);
np1p_norm = np1p./sum(np1p);
np2_norm = np2./sum(np2);
np3_norm = np3./sum(np3);

