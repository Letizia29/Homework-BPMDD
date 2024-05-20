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


%% IMPORTANT VARIABLES EXTRACTION and MANAGING MISSING VALUES
%% - Summing-up variables: summing scores for each part 

% HEALTHY CONTROLS ------------------------------------------------------
% Find missing data
np_test.np1r.HC.idx_nan = find(isnan(data.NP1RTOT(idx_HC, :)));
np_test.np1p.HC.idx_nan = find(isnan(data.NP1PTOT(idx_HC, :)));
np_test.np2.HC.idx_nan = find(isnan(data.NP2PTOT(idx_HC, :)));
np_test.np3.HC.idx_nan = find(isnan(data.NP3TOT(idx_HC, :)));
np_test.np4.HC.idx_nan = find(string(data.NP4TOT(idx_HC, :)) =='NA'); 

np_test.idx_nan.HC = union(np_test.np1r.HC.idx_nan, union(np_test.np1p.HC.idx_nan,union(np_test.np2.HC.idx_nan,np_test.np3.HC.idx_nan)));

data_without_nan_np = data(idx_HC,:);
data_without_nan_np(np_test.idx_nan.HC,:) = [];

np_test.np1r.HC.data = data_without_nan_np.NP1RTOT;
np_test.np1p.HC.data = data_without_nan_np.NP1PTOT;
np_test.np2.HC.data  = data_without_nan_np.NP2PTOT;
np_test.np3.HC.data = data_without_nan_np.NP3TOT;
np_test.np4.HC.data = data_without_nan_np.NP4TOT;

% Normalization
np_test.np1r.HC.data_normalized  = np_test.np1r.HC.data ./sum(np_test.np1r.HC.data );
np_test.np1p.HC.data_normalized = np_test.np1p.HC.data./sum(np_test.np1p.HC.data);
np_test.np2.HC.data_normalized  = np_test.np2.HC.data./sum(np_test.np2.HC.data);
np_test.np3.HC.data_normalized = np_test.np3.HC.data./sum(np_test.np3.HC.data);

%%%%%%%%%%
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
%%%%%%%%%

% PD PATIENTS ------------------------------------------------------
% Find missing data
np_test.np1r.PD.idx_nan = find(isnan(data.NP1RTOT([idx_Prodromal; idx_PD; idx_SWEDD],:)));
np_test.np1p.PD.idx_nan = find(isnan(data.NP1PTOT([idx_Prodromal; idx_PD; idx_SWEDD],:)));
np_test.np2.PD.idx_nan = find(isnan(data.NP2PTOT([idx_Prodromal; idx_PD; idx_SWEDD],:)));
np_test.np3.PD.idx_nan = find(isnan(data.NP3TOT([idx_Prodromal; idx_PD; idx_SWEDD],:)));
np_test.np4.PD.idx_nan = find(string(data.NP4TOT([idx_Prodromal; idx_PD; idx_SWEDD],:)) =='NA'); 

np_test.idx_nan.PD = union(np_test.np1r.PD.idx_nan, union(np_test.np1p.PD.idx_nan, union(np_test.np2.PD.idx_nan,np_test.np3.PD.idx_nan)));

data_without_nan_np = data([idx_Prodromal; idx_PD; idx_SWEDD],:);
data_without_nan_np(np_test.idx_nan.PD,:) = [];

np_test.np1r.PD.data = data_without_nan_np.NP1RTOT;
np_test.np1p.PD.data = data_without_nan_np.NP1PTOT;
np_test.np2.PD.data = data_without_nan_np.NP2PTOT;
np_test.np3.PD.data = data_without_nan_np.NP3TOT;
np_test.np4.PD.data = data_without_nan_np.NP4TOT;

% Normalization
np_test.np1r.PD.data_normalized = np_test.np1r.PD.data./sum(np_test.np1r.PD.data);
np_test.np1p.PD.data_normalized = np_test.np1p.PD.data./sum(np_test.np1p.PD.data);
np_test.np2.PD.data_normalized = np_test.np2.PD.data./sum(np_test.np2.PD.data);
np_test.np3.PD.data_normalized = np_test.np3.PD.data./sum(np_test.np3.PD.data);

%%%%%%%%%%
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
%%%%%%%%%%

%% - Check repetitions in PATNO - patients ID
IDs = data.PATNO;
unique_ids = unique(IDs);
if length(IDs) == length(unique_ids)
    disp('No repetitions')
end

%% - Genetic history
% genetics_hc = data_hc.GENETICS; % tutti NA
genetics_pd = data_pd.GENETICS;

familiarity_hc = data_hc.ANYFAMPD; % alcuni hanno familiarità
familiarity_pd = data_pd.ANYFAMPD; % 47% NAN 

%% - Demographics data
ethnicity_hc = data_hc.ETHNICITY; % nessun NA 

idx_nan_ethnicity = find(data_pd.ETHNICITY == "NA");
data_pd(idx_nan_ethnicity,:) = [];

ethnicity_pd = data_pd.ETHNICITY;

sex_hc = data_hc.SEX; % no NA
sex_pd = data_pd.SEX; % no NA

age_hc = data_hc.ENROLL_AGE; % già sistemati NA
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
hand_hc = data_hc.HANDED; % no NA

idx_nan_hand= find(data_pd.HANDED == "NA");
data_pd(idx_nan_hand,:) = [];

hand_pd = data_pd.HANDED; 

idx_right_hand_hc = find(hand_hc == "Right");
idx_left_hand_hc = find(hand_hc == "Left");
idx_mixed_hand_hc = find(hand_hc == "Mixed");

idx_right_hand_pd = find(hand_pd == "Right");
idx_left_hand_pd = find(hand_pd == "Left");
idx_mixed_hand_pd = find(hand_pd == "Mixed");

%% - Primary diagnosis
% SWEDD patients exclusion
ind_97_hc = find(data_hc.PRIMDIAG == 97);
ind_97_pd = find(data_pd.PRIMDIAG == 97);
for i = 1:length(ind_97_hc)
    data_hc(ind_97_hc(i),"PRIMDIAG") = {"NA"};
end
for i = 1:length(ind_97_pd)
    data_pd(ind_97_pd(i),"PRIMDIAG") = {"NA"};
end

prim_diag_hc = data_hc.PRIMDIAG;
prim_diag_pd = data_pd.PRIMDIAG;

if isempty(find(prim_diag_hc == 97))
    disp('Successful SWEDD exclusion in HC')
end
if isempty(find(prim_diag_pd == 97))
    disp('Successful SWEDD exclusion in PD')
end

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

% Right = |(right - left putamen SBR/right + left putamen SBR)| > 0.20
% Left = |(right - left putamen SBR/right + left putamen SBR)| > 0.20
% From: Ipsilateral deficits of dopaminergic neurotransmission in Parkinson s disease

% DATSCAN lateralization PD ------------------------------------------
DATSCAN_CAUDATE_lat_pd = (data_pd.DATSCAN_CAUDATE_R - data_pd.DATSCAN_CAUDATE_L)./(data_pd.DATSCAN_CAUDATE_R + data_pd.DATSCAN_CAUDATE_L);
DATSCAN_PUTAMEN_lat_pd = (data_pd.DATSCAN_PUTAMEN_R - data_pd.DATSCAN_PUTAMEN_L)./(data_pd.DATSCAN_PUTAMEN_R + data_pd.DATSCAN_PUTAMEN_L);
DATSCAN_PUTAMEN_ANT_lat_pd = (data_pd.DATSCAN_PUTAMEN_R_ANT - data_pd.DATSCAN_PUTAMEN_L_ANT)./(data_pd.DATSCAN_PUTAMEN_R_ANT + data_pd.DATSCAN_PUTAMEN_L_ANT);

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.

% With significant difference: 20%

lat_index.CAUDATE.PD.right = find(DATSCAN_CAUDATE_lat_pd>0.2);
pd_right_lat_caudate = data_pd(lat_index.CAUDATE.PD.right,:);

lat_index.CAUDATE.PD.left = find(DATSCAN_CAUDATE_lat_pd<-0.2);
pd_left_lat_caudate = data_pd(lat_index.CAUDATE.PD.left,:);

lat_index.CAUDATE.PD.none = find(DATSCAN_CAUDATE_lat_pd>-0.2 | DATSCAN_CAUDATE_lat_pd<0.2);
pd_no_lat_caudate = data_pd(lat_index.CAUDATE.PD.none,:);

lat_index.PUTAMEN.PD.right = find(DATSCAN_PUTAMEN_lat_pd>0.2);
pd_right_lat_putamen = data_pd(lat_index.PUTAMEN.PD.right,:);

lat_index.PUTAMEN.PD.left = find(DATSCAN_PUTAMEN_lat_pd<-0.2);
pd_left_lat_putamen = data_pd(lat_index.PUTAMEN.PD.left,:);

lat_index.PUTAMEN.PD.none = find(DATSCAN_PUTAMEN_lat_pd<-0.2 | DATSCAN_PUTAMEN_lat_pd<0.2);
pd_no_lat_putamen = data_pd(lat_index.PUTAMEN.PD.none,:);

lat_index.PUTAMEN_ANT.PD.right = find(DATSCAN_PUTAMEN_ANT_lat_pd>0.2);
pd_right_lat_putamen_ant = data_pd(lat_index.PUTAMEN_ANT.PD.right ,:);

lat_index.PUTAMEN_ANT.PD.left = find(DATSCAN_PUTAMEN_ANT_lat_pd<-0.2);
pd_left_lat_putamen_ant = data_pd(lat_index.PUTAMEN_ANT.PD.left,:);

lat_index.PUTAMEN_ANT.PD.none = find(DATSCAN_PUTAMEN_ANT_lat_pd<-0.2 | DATSCAN_PUTAMEN_ANT_lat_pd<0.2);
pd_no_lat_putamen_ant = data_pd(lat_index.PUTAMEN_ANT.PD.none,:);

% DATSCAN lateralization HC ------------------------------------------
DATSCAN_CAUDATE_lat_hc = (data_hc.DATSCAN_CAUDATE_R - data_hc.DATSCAN_CAUDATE_L)./(data_hc.DATSCAN_CAUDATE_R + data_hc.DATSCAN_CAUDATE_L);
DATSCAN_PUTAMEN_lat_hc = (data_hc.DATSCAN_PUTAMEN_R - data_hc.DATSCAN_PUTAMEN_L)./(data_hc.DATSCAN_PUTAMEN_R + data_hc.DATSCAN_PUTAMEN_L);
DATSCAN_PUTAMEN_ANT_lat_hc = (data_hc.DATSCAN_PUTAMEN_R_ANT - data_hc.DATSCAN_PUTAMEN_L_ANT)./(data_hc.DATSCAN_PUTAMEN_R_ANT + data_hc.DATSCAN_PUTAMEN_L_ANT);

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.

% With significant difference: 20%

lat_index.CAUDATE.HC.right = find(DATSCAN_CAUDATE_lat_hc>0.2);
hc_right_lat_caudate = data_hc(lat_index.CAUDATE.HC.right,:);

lat_index.CAUDATE.HC.left = find(DATSCAN_CAUDATE_lat_hc<-0.2);
hc_left_lat_caudate = data_hc(lat_index.CAUDATE.HC.left,:);

lat_index.CAUDATE.HC.none = find(DATSCAN_CAUDATE_lat_hc>-0.2 | DATSCAN_CAUDATE_lat_hc<0.2);
hc_no_lat_caudate = data_hc(lat_index.CAUDATE.HC.none,:);

lat_index.PUTAMEN.HC.right = find(DATSCAN_PUTAMEN_lat_hc>0.2);
hc_right_lat_putamen = data_hc(lat_index.PUTAMEN.HC.right,:);

lat_index.PUTAMEN.HC.left = find(DATSCAN_PUTAMEN_lat_hc<-0.2);
hc_left_lat_putamen = data_hc(lat_index.PUTAMEN.HC.left,:);

lat_index.PUTAMEN.HC.none = find(DATSCAN_PUTAMEN_lat_hc<-0.2 | DATSCAN_PUTAMEN_lat_hc<0.2);
hc_no_lat_putamen = data_hc(lat_index.PUTAMEN.HC.none,:);

lat_index.PUTAMEN_ANT.HC.right = find(DATSCAN_PUTAMEN_ANT_lat_hc>0.2);
hc_right_lat_putamen_ant = data_hc(lat_index.PUTAMEN_ANT.HC.right ,:);

lat_index.PUTAMEN_ANT.HC.left = find(DATSCAN_PUTAMEN_ANT_lat_hc<-0.2);
hc_left_lat_putamen_ant = data_hc(lat_index.PUTAMEN_ANT.HC.left,:);

lat_index.PUTAMEN_ANT.HC.none = find(DATSCAN_PUTAMEN_ANT_lat_hc<-0.2 | DATSCAN_PUTAMEN_ANT_lat_hc<0.2);
hc_no_lat_putamen_ant = data_hc(lat_index.PUTAMEN_ANT.HC.none,:);


%%%%%%%%%%%%%% Graphical comparison
% graphical measures
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
