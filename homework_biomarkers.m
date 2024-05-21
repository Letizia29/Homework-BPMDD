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

idx_nan.ethnicity = find(data_pd.ETHNICITY == "NA");
data_pd(idx_nan.ethnicity,:) = [];

ethnicity_pd = data_pd.ETHNICITY;

sex_hc = data_hc.SEX; % no NA
sex_pd = data_pd.SEX; % no NA

age_hc = data_hc.ENROLL_AGE; % già sistemati NA
age_pd = data_pd.ENROLL_AGE;

% computing missing age data
idx_nan.age.HC = find(isnan(age_hc));
idx_nan.age.PD = find(isnan(age_pd));

% age at the time of the DAT SCAN ---- HC
for i = 1:length(idx_nan.age.HC)
    idx_nan_temp = idx_nan.age.HC(i);
    datscan_date = char(data_hc.DATSCAN_DATE(idx_nan_temp));
    birth_date = char(data_hc.BIRTHDT(idx_nan_temp));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    age_hc(idx_nan_temp) = datscan_year - birth_year;
end

% age at the time of the DAT SCAN ---- PD
for i = 1:length(idx_nan.age.PD)
    idx_nan_temp = idx_nan.age.PD(i);
    datscan_date = char(data_pd.DATSCAN_DATE(idx_nan_temp));
    birth_date = char(data_pd.BIRTHDT(idx_nan_temp));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    age_pd(idx_nan_temp) = datscan_year - birth_year;
end

%% - Dominant hand
hand_hc = data_hc.HANDED; % no NA

idx_nan.hand= find(data_pd.HANDED == "NA");
data_pd(idx_nan.hand,:) = [];

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

if isempty(find(prim_diag_hc == 97,1))
    disp('Successful SWEDD exclusion in HC')
end
if isempty(find(prim_diag_pd == 97,1))
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

% Remove not completed DAT SCANS ( DATSCAN == 0 )
idx_datscan_not_compl_pd = find(data_pd.DATSCAN == 0);
data_pd(idx_datscan_not_compl_pd,:) = [];

idx_datscan_not_compl_hc = find(data_hc.DATSCAN == 0);
data_hc(idx_datscan_not_compl_hc,:) = [];


% Right = |(right - left putamen SBR/right + left putamen SBR)| > 0.20
% Left = |(right - left putamen SBR/right + left putamen SBR)| > 0.20
% From: Ipsilateral deficits of dopaminergic neurotransmission in Parkinson s disease

% DATSCAN lateralization PD ------------------------------------------
DATSCAN.CAUDATE_lat.PD = (data_pd.DATSCAN_CAUDATE_R - data_pd.DATSCAN_CAUDATE_L)./(data_pd.DATSCAN_CAUDATE_R + data_pd.DATSCAN_CAUDATE_L);
DATSCAN.PUTAMEN_lat.PD = (data_pd.DATSCAN_PUTAMEN_R - data_pd.DATSCAN_PUTAMEN_L)./(data_pd.DATSCAN_PUTAMEN_R + data_pd.DATSCAN_PUTAMEN_L);
DATSCAN.PUTAMEN_ANT_lat.PD = (data_pd.DATSCAN_PUTAMEN_R_ANT - data_pd.DATSCAN_PUTAMEN_L_ANT)./(data_pd.DATSCAN_PUTAMEN_R_ANT + data_pd.DATSCAN_PUTAMEN_L_ANT);

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.

% With significant difference: 20%

lateralization.CAUDATE.PD.right.index = find(DATSCAN.CAUDATE_lat.PD>0.2);
lateralization.CAUDATE.PD.right.population = data_pd(lateralization.CAUDATE.PD.right.index,:);

lateralization.CAUDATE.PD.left.index = find(DATSCAN.CAUDATE_lat.PD<-0.2);
lateralization.CAUDATE.PD.left.population = data_pd(lateralization.CAUDATE.PD.left.index,:);

lateralization.CAUDATE.PD.none.index = find(DATSCAN.CAUDATE_lat.PD>-0.2 | DATSCAN.CAUDATE_lat.PD<0.2);
lateralization.CAUDATE.PD.none.population = data_pd(lateralization.CAUDATE.PD.none.index,:);

lateralization.PUTAMEN.PD.right.index = find(DATSCAN.PUTAMEN_lat.PD>0.2);
lateralization.PUTAMEN.PD.right.population = data_pd(lateralization.PUTAMEN.PD.right.index,:);

lateralization.PUTAMEN.PD.left.index = find(DATSCAN.PUTAMEN_lat.PD<-0.2);
lateralization.PUTAMEN.PD.left.population = data_pd(lateralization.PUTAMEN.PD.left.index,:);

lateralization.PUTAMEN.PD.none.index = find(DATSCAN.PUTAMEN_lat.PD<-0.2 | DATSCAN.PUTAMEN_lat.PD<0.2);
lateralization.PUTAMEN.PD.none.population = data_pd(lateralization.PUTAMEN.PD.none.index,:);

lateralization.PUTAMEN_ANT.PD.right.index = find(DATSCAN.PUTAMEN_ANT_lat.PD>0.2);
lateralization.PUTAMEN_ANT.PD.right.population = data_pd(lateralization.PUTAMEN_ANT.PD.right.index ,:);

lateralization.PUTAMEN_ANT.PD.left.index = find(DATSCAN.PUTAMEN_ANT_lat.PD<-0.2);
lateralization.PUTAMEN_ANT.PD.left.population = data_pd(lateralization.PUTAMEN_ANT.PD.left.index,:);

lateralization.PUTAMEN_ANT.PD.none.index = find(DATSCAN.PUTAMEN_ANT_lat.PD<-0.2 | DATSCAN.PUTAMEN_ANT_lat.PD<0.2);
lateralization.PUTAMEN_ANT.PD.none.population = data_pd(lateralization.PUTAMEN_ANT.PD.none.index,:);

% DATSCAN lateralization HC ------------------------------------------
DATSCAN.CAUDATE_lat.HC = (data_hc.DATSCAN_CAUDATE_R - data_hc.DATSCAN_CAUDATE_L)./(data_hc.DATSCAN_CAUDATE_R + data_hc.DATSCAN_CAUDATE_L);
DATSCAN.PUTAMEN_lat.HC = (data_hc.DATSCAN_PUTAMEN_R - data_hc.DATSCAN_PUTAMEN_L)./(data_hc.DATSCAN_PUTAMEN_R + data_hc.DATSCAN_PUTAMEN_L);
DATSCAN.PUTAMEN_ANT_lat.HC = (data_hc.DATSCAN_PUTAMEN_R_ANT - data_hc.DATSCAN_PUTAMEN_L_ANT)./(data_hc.DATSCAN_PUTAMEN_R_ANT + data_hc.DATSCAN_PUTAMEN_L_ANT);

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.

% With significant difference: 20%

lateralization.CAUDATE.HC.right.index = find(DATSCAN.CAUDATE_lat.HC>0.2);
lateralization.CAUDATE.HC.right.population = data_hc(lateralization.CAUDATE.HC.right.index,:);

lateralization.CAUDATE.HC.left.index = find(DATSCAN.CAUDATE_lat.HC<-0.2);
lateralization.CAUDATE.HC.left.population = data_hc(lateralization.CAUDATE.HC.left.index,:);

lateralization.CAUDATE.HC.none.index = find(DATSCAN.CAUDATE_lat.HC>-0.2 | DATSCAN.CAUDATE_lat.HC<0.2);
lateralization.CAUDATE.HC.none.population = data_hc(lateralization.CAUDATE.HC.none.index,:);

lateralization.PUTAMEN.HC.right.index = find(DATSCAN.PUTAMEN_lat.HC>0.2);
lateralization.PUTAMEN.HC.right.population = data_hc(lateralization.PUTAMEN.HC.right.index,:);

lateralization.PUTAMEN.HC.left.index = find(DATSCAN.PUTAMEN_lat.HC<-0.2);
lateralization.PUTAMEN.HC.left.population = data_hc(lateralization.PUTAMEN.HC.left.index,:);

lateralization.PUTAMEN.HC.none.index = find(DATSCAN.PUTAMEN_lat.HC<-0.2 | DATSCAN.PUTAMEN_lat.HC<0.2);
lateralization.PUTAMEN.HC.none.population = data_hc(lateralization.PUTAMEN.HC.none.index,:);

lateralization.PUTAMEN_ANT.HC.right.index = find(DATSCAN.PUTAMEN_ANT_lat.HC>0.2);
lateralization.PUTAMEN_ANT.HC.right.population = data_hc(lateralization.PUTAMEN_ANT.HC.right.index ,:);

lateralization.PUTAMEN_ANT.HC.left.index = find(DATSCAN.PUTAMEN_ANT_lat.HC<-0.2);
lateralization.PUTAMEN_ANT.HC.left.population = data_hc(lateralization.PUTAMEN_ANT.HC.left.index,:);

lateralization.PUTAMEN_ANT.HC.none.index = find(DATSCAN.PUTAMEN_ANT_lat.HC<-0.2 | DATSCAN.PUTAMEN_ANT_lat.HC<0.2);
lateralization.PUTAMEN_ANT.HC.none.population = data_hc(lateralization.PUTAMEN_ANT.HC.none.index,:);


%%%%%%%%%%%%%% Graphical comparison
% graphical measures
% x_left = [0; 0; 1556; 1556];
% y_left = [0; 1.9600; 1.9600; 0];
% x_right = [0; 0; 1556; 1556];
% y_right = [-1.9600; 0; 0; -1.9600];
% 
% figure(4), hold on, axis tight
% plot(DATSCAN_CAUDATE_lat.PD, 'k')
% title('Caudate DAT signal lateralization (right-left)')
% ylabel('[adim]')
% xlabel('Subjects')
% yline(0, 'r', 'LineWidth', 2)
% patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
% patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
% legend('DAT_{right} - DAT_{left}', 'Threshold', 'DAT_{right} bigger', 'DAT_{left} bigger')
% 
% figure(5), hold on, axis tight
% plot(DATSCAN_PUTAMEN_lat.PD, 'k')
% title('Putamen DAT signal lateralization (right-left)')
% ylabel('[adim]')
% xlabel('Subjects')
% yline(0, 'r', 'LineWidth', 2)
% patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
% patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
% legend('DAT_{right} - DAT_{left}', 'Threshold', 'DAT_{right} bigger', 'DAT_{left} bigger')
% 
% figure(6), hold on, axis tight
% plot(DATSCAN_PUTAMEN_ANT_lat.PD, 'k')
% title('Putamen ant. DAT signal lateralization (right-left)')
% ylabel('[adim]')
% xlabel('Subjects')
% yline(0, 'r', 'LineWidth', 2)
% patch(x_left, y_left, 'b', 'FaceAlpha', 0.3)
% patch(x_right, y_right, 'r', 'FaceAlpha', 0.3)
% legend('DAT_{right} - DAT_{left}', 'Threshold', 'DAT_{right} bigger', 'DAT_{left} bigger')
%%%%%%%%%%%%%%%%%

%% PATTERN MISSING VALUES

% np1r, np1p, np2, np3, np4, genetics, familiarity, ethnicity, sex, age,
% hand, prim_diag

missing_values.pd = zeros(size(data_pd, 1), 12);
missing_values.hc = zeros(size(data_hc, 1), 12);

missing_values.pd(np_test.np1p.PD.idx_nan, 1) = 1;
missing_values.hc(np_test.np1p.HC.idx_nan, 1) = 1;

missing_values.pd(np_test.np1r.PD.idx_nan, 2) = 1;
missing_values.hc(np_test.np1r.HC.idx_nan, 2) = 1;

missing_values.pd(np_test.np2.PD.idx_nan, 3) = 1;
missing_values.hc(np_test.np2.HC.idx_nan, 3) = 1;

missing_values.pd(np_test.np3.PD.idx_nan, 4) = 1;
missing_values.hc(np_test.np3.HC.idx_nan, 4) = 1;

missing_values.pd(np_test.np4.PD.idx_nan, 5) = 1;
missing_values.hc(np_test.np4.HC.idx_nan, 5) = 1;

missing_values.pd(find(data_pd.GENETICS == "NA"), 6) = 1;
missing_values.hc(find(data_hc.GENETICS == "NA"), 6) = 1;

missing_values.pd(find(data_pd.ANYFAMPD == "NA"), 7) = 1;
missing_values.hc(find(data_hc.ANYFAMPD == "NA"), 7) = 1;

missing_values.pd(find(data_pd.ETHNICITY == "NA"), 8) = 1;
missing_values.hc(find(data_hc.ETHNICITY == "NA"), 8) = 1;

missing_values.pd(find(data_pd.SEX == "NA"), 9) = 1;
missing_values.hc(find(data_hc.SEX == "NA"), 9) = 1;

missing_values.pd(find(isnan(data_pd.ENROLL_AGE)), 10) = 1;
missing_values.hc(find(isnan(data_hc.ENROLL_AGE)), 10) = 1;

missing_values.pd(find(data_pd.HANDED == "NA"), 11) = 1;
missing_values.hc(find(data_hc.HANDED == "NA"), 11) = 1;

missing_values.pd(find(string(data_pd.PRIMDIAG) == "NA"), 12) = 1;
missing_values.hc(find(string(data_hc.PRIMDIAG) == "NA"), 12) = 1;


figure
subplot(121)
imagesc(missing_values.hc)
cmap = jet(2);
colormap(cmap)
hold on
L = line(ones(2), ones(2), 'LineWidth', 2);
set(L, {'color'}, mat2cell(cmap, ones(1, 2), 3));
title('HC')
xticks(1:12)
xticklabels({'np1r', 'np1p', 'np2', 'np3', 'np4', 'genetics', 'familiarity', 'ethnicity', 'sex', 'age', 'hand', 'prim_diag'})
legend('Data', 'NaN')


subplot(122)
imagesc(missing_values.pd)
title('PD')
cmap = jet(2);
colormap(cmap)
hold on
L = line(ones(2), ones(2), 'LineWidth', 2);
set(L, {'color'}, mat2cell(cmap, ones(1, 2), 3));
xticks(1:12)
xticklabels({'np1r', 'np1p', 'np2', 'np3', 'np4', 'genetics', 'familiarity', 'ethnicity', 'sex', 'age', 'hand', 'prim_diag'})
legend('Data', 'NaN')


