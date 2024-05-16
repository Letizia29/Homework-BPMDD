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

% division HC - PD
data_pd = data([idx_Prodromal; idx_PD; idx_SWEDD],:);
data_hc = data(idx_HC, :);

%%%%%%%%%%% histograms
% figure(1)
% for i=1:6
%     subplot(2,3,i)
%     histogram(NON_HC_DAT_SCAN_SPECT(:,i),'FaceColor','auto','Normalization','probability')
%     xlabel('Striatal binding ratio [adim]')
%     title(['SBR in ' ROIs_labels(i)])
%     hold on 
%     mu_non_hc = mean(NON_HC_DAT_SCAN_SPECT(:,i));
%     xline(mu_non_hc,'LineWidth',2,'Color','b')
%     histogram(HC_DAT_SCAN_SPECT(:,i),'FaceColor','auto','Normalization','probability')
%     mu_hc = mean(HC_DAT_SCAN_SPECT(:,i));
%     xline(mu_hc,'LineWidth',2,'Color','r')
%     xlabel('Striatal binding ratio [adim]')
%     title(['SBR in ' ROIs_labels(i)])
%     hold off
%     legend('Non healthy','Average DAT non healthy','Healthy controls','Average DAT healthy')
% end
%%%%%%%%%%%%


%% IMPORTANT VARIABLES EXTRACTION:
%% - Summing-up variables: summing scores for each part 

% HEALTHY CONTROLS ------------------------------------------------------
np1r_HC = data.NP1RTOT(idx_HC, :);
np1p_HC = data.NP1PTOT(idx_HC, :);
np2_HC = data.NP2PTOT(idx_HC, :);
np3_HC = data.NP3TOT(idx_HC, :);
np4_HC = data.NP4TOT(idx_HC, :);

% Find missing data
idx_nan_np1r_HC = find(isnan(np1r_HC));
idx_nan_np1p_HC = find(isnan(np1p_HC));
idx_nan_np2_HC = find(isnan(np2_HC));
idx_nan_np3_HC = find(isnan(np3_HC));
idx_nan_np4_HC = find(string(np4_HC) =='NA'); 
idx_nan_np1_HC = union(idx_nan_np1r_HC, idx_nan_np1p_HC);
idx_nan_np23_HC = union(idx_nan_np2_HC, idx_nan_np3_HC);
idx_nan_np123_HC = union(idx_nan_np23_HC, idx_nan_np1_HC);
idx_nan_np_HC = union(idx_nan_np123_HC, idx_nan_np4_HC);

data_without_nan_np = data(idx_HC,:);
data_without_nan_np(idx_nan_np123_HC,:) = [];

np1r_HC = data_without_nan_np.NP1RTOT;
np1p_HC = data_without_nan_np.NP1PTOT;
np2_HC = data_without_nan_np.NP2PTOT;
np3_HC = data_without_nan_np.NP3TOT;
np4_HC = data_without_nan_np.NP4TOT;

% Normalization
np1r_norm_HC = np1r_HC./sum(np1r_HC);
np1p_norm_HC = np1p_HC./sum(np1p_HC);
np2_norm_HC = np2_HC./sum(np2_HC);
np3_norm_HC = np3_HC./sum(np3_HC);

% figure(2), hold on
% histogram(np1r_norm_HC, "FaceColor", 'r')
% histogram(np1p_norm_HC, "FaceColor", 'b')
% histogram(np2_norm_HC, "FaceColor", 'g')
% histogram(np3_norm_HC, "FaceColor", 'm')
% axis tight
% legend('NP1R', 'NP1P', 'NP2', 'NP3')
% title('Normalized sum of symptoms scores - HC')
% xlabel('Normalized score')
% ylabel('Subjects')

% PD PATIENTS ------------------------------------------------------
np1r_PD = data.NP1RTOT([idx_Prodromal; idx_PD; idx_SWEDD],:);
np1p_PD = data.NP1PTOT([idx_Prodromal; idx_PD; idx_SWEDD],:);
np2_PD = data.NP2PTOT([idx_Prodromal; idx_PD; idx_SWEDD],:);
np3_PD = data.NP3TOT([idx_Prodromal; idx_PD; idx_SWEDD],:);
np4_PD = data.NP4TOT([idx_Prodromal; idx_PD; idx_SWEDD],:);

% Find missing data
idx_nan_np1r_PD = find(isnan(np1r_PD));
idx_nan_np1p_PD = find(isnan(np1p_PD));
idx_nan_np2_PD = find(isnan(np2_PD));
idx_nan_np3_PD = find(isnan(np3_PD));
idx_nan_np4_PD = find(string(np4_PD) =='NA'); 
idx_nan_np1_PD = union(idx_nan_np1r_PD, idx_nan_np1p_PD);
idx_nan_np23_PD = union(idx_nan_np2_PD, idx_nan_np3_PD);
idx_nan_np123_PD = union(idx_nan_np23_PD, idx_nan_np1_PD);
idx_nan_np_PD = union(idx_nan_np123_PD, idx_nan_np4_PD);

data_without_nan_np = data([idx_Prodromal; idx_PD; idx_SWEDD],:);
data_without_nan_np(idx_nan_np123_PD,:) = [];

np1r_PD = data_without_nan_np.NP1RTOT;
np1p_PD = data_without_nan_np.NP1PTOT;
np2_PD = data_without_nan_np.NP2PTOT;
np3_PD = data_without_nan_np.NP3TOT;
np4_PD = data_without_nan_np.NP4TOT;

% Normalization
np1r_norm_PD = np1r_PD./sum(np1r_PD);
np1p_norm_PD = np1p_PD./sum(np1p_PD);
np2_norm_PD = np2_PD./sum(np2_PD);
np3_norm_PD = np3_PD./sum(np3_PD);

% figure(3), hold on
% histogram(np1r_norm_PD, "FaceColor", 'r')
% histogram(np1p_norm_PD, "FaceColor", 'b')
% histogram(np2_norm_PD, "FaceColor", 'g')
% histogram(np3_norm_PD, "FaceColor", 'm')
% axis tight
% legend('NP1R', 'NP1P', 'NP2', 'NP3')
% title('Normalized sum of symptoms scores - PD')
% xlabel('Normalized score')
% ylabel('Subjects')

%% - Check repetitions in PATNO - patients ID
IDs = data.PATNO;
unique_ids = unique(IDs);
if length(IDs) == length(unique_ids)
    disp('No repetitions')
end

%% - Genetic history
genetics_hc = data_hc.GENETICS;
genetics_pd = data_pd.GENETICS;

familiarity_hc = data_hc.ANYFAMPD;
familiarity_pd = data_pd.ANYFAMPD;

%% - Demographics data
ethnicity_hc = data_hc.ETHNICITY;
ethnicity_pd = data_pd.ETHNICITY;

sex_hc = data_hc.SEX;
sex_pd = data_pd.SEX;

age_hc = data_hc.ENROLL_AGE;
age_pd = data_pd.ENROLL_AGE;

% computing missing age data
idx_nan_age_hc = find(isnan(age_hc));
idx_nan_age_pd = find(isnan(age_pd));

% age at the time of the DAT SCAN ---- HC
for i = 1:length(idx_nan_age_hc)
    idx_nan = idx_nan_age_hc(i);
    datscan_date = char(data_hc.DATSCAN_DATE(idx_nan));
    birth_date = char(data_hc.BIRTHDT(idx_nan));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    age_hc(idx_nan) = datscan_year - birth_year;
end

% age at the time of the DAT SCAN ---- PD
for i = 1:length(idx_nan_age_pd)
    idx_nan = idx_nan_age_pd(i);
    datscan_date = char(data_pd.DATSCAN_DATE(idx_nan));
    birth_date = char(data_pd.BIRTHDT(idx_nan));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    age_pd(idx_nan) = datscan_year - birth_year;
end

%% - Dominant hand
hand_hc = data_hc.HANDED;
hand_pd = data_pd.HANDED;

idx_right_hand_hc = find(hand_hc == "Right");
idx_left_hand_hc = find(hand_hc == "Left");
idx_mixed_hand_hc = find(hand_hc == "Mixed");

idx_right_hand_pd = find(hand_pd == "Right");
idx_left_hand_pd = find(hand_pd == "Left");
idx_mixed_hand_pd = find(hand_pd == "Mixed");

%% - Primary diagnosis
prim_diag_hc = data_hc.PRIMDIAG;
prim_diag_pd = data_pd.PRIMDIAG;



%% MISSING VALUES MANAGING
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

%% DAT SCAN lateralization
DATSCAN_CAUDATE_lat = data.DATSCAN_CAUDATE_R - data.DATSCAN_CAUDATE_L;
DATSCAN_PUTAMEN_lat = data.DATSCAN_PUTAMEN_R - data.DATSCAN_PUTAMEN_L;
DATSCAN_PUTAMEN_ANT_lat = data.DATSCAN_PUTAMEN_R_ANT - data.DATSCAN_PUTAMEN_L_ANT;

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.

%%%%%%%%%%%%%% Graphical comparison
% % graphical measures
% x_left = [0; 0; 1556; 1556];
% y_left = [0; 1.9600; 1.9600; 0];
% x_right = [0; 0; 1556; 1556];
% y_right = [-1.9600; 0; 0; -1.9600];
% 
% figure(4), hold on, axis tight
% plot(DATSCAN_CAUDATE_lat, 'k')
% title('Caudate DAT signal lateralization (right-left)')
% ylabel('[adim]')
% xlabel('Subjects')
% yline(0, 'r', 'LineWidth', 2)
% patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
% patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
% legend('DAT_{right} - DAT_{left}', 'Threshold', 'DAT_{right} bigger', 'DAT_{left} bigger')
% 
% figure(5), hold on, axis tight
% plot(DATSCAN_PUTAMEN_lat, 'k')
% title('Putamen DAT signal lateralization (right-left)')
% ylabel('[adim]')
% xlabel('Subjects')
% yline(0, 'r', 'LineWidth', 2)
% patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
% patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
% legend('DAT_{right} - DAT_{left}', 'Threshold', 'DAT_{right} bigger', 'DAT_{left} bigger')
% 
% figure(6), hold on, axis tight
% plot(DATSCAN_PUTAMEN_ANT_lat, 'k')
% title('Putamen ant. DAT signal lateralization (right-left)')
% ylabel('[adim]')
% xlabel('Subjects')
% yline(0, 'r', 'LineWidth', 2)
% patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
% patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
% legend('DAT_{right} - DAT_{left}', 'Threshold', 'DAT_{right} bigger', 'DAT_{left} bigger')
%%%%%%%%%%%%%%%%%
