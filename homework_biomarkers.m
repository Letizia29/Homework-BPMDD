clear all
close all
clc

%% load data
data = readtable("Patient_Master.csv");

DATSCAN.TOT = [data.DATSCAN_CAUDATE_R data.DATSCAN_CAUDATE_L data.DATSCAN_PUTAMEN_R  data.DATSCAN_PUTAMEN_L data.DATSCAN_PUTAMEN_R_ANT data.DATSCAN_PUTAMEN_L_ANT];
%DAT_SCAN_PET =  [data.AV133_RCAUD_S data.AV133_LCAUD_S data.AV133_RPUTANT_S data.AV133_RPUTPOST_S data.AV133_LPUTANT_S data.AV133_LPUTPOST_S];

ROIs_labels = ["Right Caudate", "Left Caudate", "Right Putamen", "Left Putamen", "Left Anterior Putamen", "Right Anterior Putamen"];
idx_variables.HC = find(string(data.COHORT)=='HC');
idx_variables.PD = find(string(data.COHORT)=='PD');
idx_variables.SWEDD = find(string(data.COHORT)=='SWEDD');
idx_variables.Prodromal = find(string(data.COHORT)=='Prodromal');


% division HC - PD
data_pd = data([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:);
data_hc = data(idx_variables.HC, :);

%% SAVE NEW DATASET WITH NANs for plotting in R
new_data_hc = table();
new_data_pd = table();

important_var = [34, 41, 55, 94, 104, 5, 11, 7, 8, 4, 153, 154, 9, 12];
names = data_hc.Properties.VariableNames;
% np1ptot, np1rtot, np2tot, np3tot, np4tot, genetics, familiarity, ethnicity, sex, age,
% height, weight, hand, primary diagnosis, 

for i = 1:length(important_var)
    new_data_hc(:,i) = data_hc(:,names(important_var(i)));
    new_data_pd(:,i) = data_pd(:,names(important_var(i)));
end

new_data_hc.Properties.VariableNames = names(1, important_var);
new_data_pd.Properties.VariableNames = names(1, important_var);

writetable(new_data_hc, 'new_data_hc.csv');
writetable(new_data_pd, 'new_data_pd.csv');

clear new_data_pd new_data_hc

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
np_test.np1r.HC.idx_nan = find(isnan(data.NP1RTOT(idx_variables.HC, :)));
np_test.np1p.HC.idx_nan = find(isnan(data.NP1PTOT(idx_variables.HC, :)));
np_test.np2.HC.idx_nan = find(isnan(data.NP2PTOT(idx_variables.HC, :)));
np_test.np3.HC.idx_nan = find(isnan(data.NP3TOT(idx_variables.HC, :)));
np_test.np4.HC.idx_nan = find(string(data.NP4TOT(idx_variables.HC, :)) =='NA'); 

np_test.idx_nan.HC = union(np_test.np1r.HC.idx_nan, union(np_test.np1p.HC.idx_nan,union(np_test.np2.HC.idx_nan,np_test.np3.HC.idx_nan)));

data_without_nan_np = data(idx_variables.HC,:);
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
np_test.np1r.PD.idx_nan = find(isnan(data.NP1RTOT([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:)));
np_test.np1p.PD.idx_nan = find(isnan(data.NP1PTOT([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:)));
np_test.np2.PD.idx_nan = find(isnan(data.NP2PTOT([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:)));
np_test.np3.PD.idx_nan = find(isnan(data.NP3TOT([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:)));
np_test.np4.PD.idx_nan = find(string(data.NP4TOT([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:)) =='NA'); 

np_test.idx_nan.PD = union(np_test.np1r.PD.idx_nan, union(np_test.np1p.PD.idx_nan, union(np_test.np2.PD.idx_nan,np_test.np3.PD.idx_nan)));

data_without_nan_np = data([idx_variables.Prodromal; idx_variables.PD; idx_variables.SWEDD],:);
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

% save in variables
variables.np_test = np_test;

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

%% Remove patients with missing symptoms data
data_hc(np_test.idx_nan.HC, :) = [];
data_pd(np_test.idx_nan.PD, :) = [];




%% - Check repetitions in PATNO - patients ID
IDs = data.PATNO;
variables.unique_ids = unique(IDs);
if length(IDs) == length(variables.unique_ids)
    disp('No repetitions')
end

%% - Genetic history
% genetics_hc = data_hc.GENETICS; % tutti NA
variables.genetics.pd = data_pd.GENETICS;

variables.familiarity.hc = data_hc.ANYFAMPD; % alcuni hanno familiarità
variables.familiarity.pd = data_pd.ANYFAMPD; % 47% NAN 

%% - Demographics data
variables.ethnicity.hc = data_hc.ETHNICITY; % nessun NA 

% idx_nan.ethnicity = find(data_pd.ETHNICITY == "NA");
% data_pd(idx_nan.ethnicity,:) = [];

variables.ethnicity.pd = data_pd.ETHNICITY;

variables.sex.hc = data_hc.SEX; % no NA
variables.sex.pd = data_pd.SEX; % no NA

variables.age.hc = data_hc.ENROLL_AGE; % già sistemati NA
variables.age.pd = data_pd.ENROLL_AGE;

% computing missing age data
idx_nan.age.HC = find(isnan(variables.age.hc));
idx_nan.age.PD = find(isnan(variables.age.pd));

% age at the time of the DAT SCAN ---- HC
for i = 1:length(idx_nan.age.HC)
    idx_nan_temp = idx_nan.age.HC(i);
    datscan_date = char(data_hc.DATSCAN_DATE(idx_nan_temp));
    birth_date = char(data_hc.BIRTHDT(idx_nan_temp));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    variables.age.hc(idx_nan_temp) = datscan_year - birth_year;
end

% age at the time of the DAT SCAN ---- PD
for i = 1:length(idx_nan.age.PD)
    idx_nan_temp = idx_nan.age.PD(i);
    datscan_date = char(data_pd.DATSCAN_DATE(idx_nan_temp));
    birth_date = char(data_pd.BIRTHDT(idx_nan_temp));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    variables.age.pd(idx_nan_temp) = datscan_year - birth_year;
end

% height
variables.height.hc = data_hc.HTCM;
variables.height.pd = data_pd.HTCM;

% weight
variables.weight.hc = data_hc.WGTKG;
variables.weight.pd = data_pd.WGTKG;


%% - Dominant hand
variables.hand.hc = data_hc.HANDED; % no NA

% idx_nan.hand= find(data_pd.HANDED == "NA");
% data_pd(idx_nan.hand,:) = [];

variables.hand.pd = data_pd.HANDED; 

idx_variables.hand.right.hc = find(variables.hand.hc == "Right");
idx_variables.hand.left.hc  = find(variables.hand.hc == "Left");
idx_variables.hand.mixed.hc  = find(variables.hand.hc == "Mixed");

idx_variables.hand.right.pd = find(variables.hand.pd == "Right");
idx_variables.hand.left.pd = find(variables.hand.pd == "Left");
idx_variables.hand.mixed.pd = find(variables.hand.pd == "Mixed");

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

variables.prim_diag.hc = data_hc.PRIMDIAG;
variables.prim_diag.pd = data_pd.PRIMDIAG;

if isempty(find(variables.prim_diag.hc == 97,1))
    disp('Successful SWEDD exclusion in HC')
end
if isempty(find(variables.prim_diag.pd == 97,1))
    disp('Successful SWEDD exclusion in PD')
end



%% MISSING VALUES MANAGING
%% Analysis missing data ex family history
answers_fam_pd = data.ANYFAMPD;
idx_variables.familiarity.yes.pd = find(string(answers_fam_pd) =='1');
%yes_fam_pd = answers_fam_pd
num_yes_fam_pd = length(idx_variables.familiarity.yes.pd);
idx_variables.familiarity.no.pd = find(string(answers_fam_pd) =='0');
num_no_fam_pd = length(idx_variables.familiarity.no.pd);
idx_nan.familiarity.pd = find(string(answers_fam_pd) =='NA');
num_nan_fam_pd = length(idx_nan.familiarity.pd);

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

lateralization.CAUDATE.PD.none.index = find(DATSCAN.CAUDATE_lat.PD>=-0.2 & DATSCAN.CAUDATE_lat.PD<=0.2);
lateralization.CAUDATE.PD.none.population = data_pd(lateralization.CAUDATE.PD.none.index,:);

lateralization.PUTAMEN.PD.right.index = find(DATSCAN.PUTAMEN_lat.PD>0.2);
lateralization.PUTAMEN.PD.right.population = data_pd(lateralization.PUTAMEN.PD.right.index,:);

lateralization.PUTAMEN.PD.left.index = find(DATSCAN.PUTAMEN_lat.PD<-0.2);
lateralization.PUTAMEN.PD.left.population = data_pd(lateralization.PUTAMEN.PD.left.index,:);

lateralization.PUTAMEN.PD.none.index = find(DATSCAN.PUTAMEN_lat.PD>=-0.2 & DATSCAN.PUTAMEN_lat.PD<=0.2);
lateralization.PUTAMEN.PD.none.population = data_pd(lateralization.PUTAMEN.PD.none.index,:);

lateralization.PUTAMEN_ANT.PD.right.index = find(DATSCAN.PUTAMEN_ANT_lat.PD>0.2);
lateralization.PUTAMEN_ANT.PD.right.population = data_pd(lateralization.PUTAMEN_ANT.PD.right.index ,:);

lateralization.PUTAMEN_ANT.PD.left.index = find(DATSCAN.PUTAMEN_ANT_lat.PD<-0.2);
lateralization.PUTAMEN_ANT.PD.left.population = data_pd(lateralization.PUTAMEN_ANT.PD.left.index,:);

lateralization.PUTAMEN_ANT.PD.none.index = find(DATSCAN.PUTAMEN_ANT_lat.PD>=-0.2 & DATSCAN.PUTAMEN_ANT_lat.PD<=0.2);
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

lateralization.CAUDATE.HC.none.index = find(DATSCAN.CAUDATE_lat.HC>=-0.2 & DATSCAN.CAUDATE_lat.HC<=0.2);
lateralization.CAUDATE.HC.none.population = data_hc(lateralization.CAUDATE.HC.none.index,:);

lateralization.PUTAMEN.HC.right.index = find(DATSCAN.PUTAMEN_lat.HC>0.2);
lateralization.PUTAMEN.HC.right.population = data_hc(lateralization.PUTAMEN.HC.right.index,:);

lateralization.PUTAMEN.HC.left.index = find(DATSCAN.PUTAMEN_lat.HC<-0.2);
lateralization.PUTAMEN.HC.left.population = data_hc(lateralization.PUTAMEN.HC.left.index,:);

lateralization.PUTAMEN.HC.none.index = find(DATSCAN.PUTAMEN_lat.HC>=-0.2 & DATSCAN.PUTAMEN_lat.HC<=0.2);
lateralization.PUTAMEN.HC.none.population = data_hc(lateralization.PUTAMEN.HC.none.index,:);

lateralization.PUTAMEN_ANT.HC.right.index = find(DATSCAN.PUTAMEN_ANT_lat.HC>0.2);
lateralization.PUTAMEN_ANT.HC.right.population = data_hc(lateralization.PUTAMEN_ANT.HC.right.index ,:);

lateralization.PUTAMEN_ANT.HC.left.index = find(DATSCAN.PUTAMEN_ANT_lat.HC<-0.2);
lateralization.PUTAMEN_ANT.HC.left.population = data_hc(lateralization.PUTAMEN_ANT.HC.left.index,:);

lateralization.PUTAMEN_ANT.HC.none.index = find(DATSCAN.PUTAMEN_ANT_lat.HC>=-0.2 & DATSCAN.PUTAMEN_ANT_lat.HC<=0.2);
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
% height, weight, hand, primary diagnosis, 
% 
% missing_values.pd = zeros(size(data_pd, 1), 14);
% missing_values.hc = zeros(size(data_hc, 1), 14);
% 
% missing_values.pd(np_test.np1p.PD.idx_nan, 1) = 1;
% missing_values.hc(np_test.np1p.HC.idx_nan, 1) = 1;
% 
% missing_values.pd(np_test.np1r.PD.idx_nan, 2) = 1;
% missing_values.hc(np_test.np1r.HC.idx_nan, 2) = 1;
% 
% missing_values.pd(np_test.np2.PD.idx_nan, 3) = 1;
% missing_values.hc(np_test.np2.HC.idx_nan, 3) = 1;
% 
% missing_values.pd(np_test.np3.PD.idx_nan, 4) = 1;
% missing_values.hc(np_test.np3.HC.idx_nan, 4) = 1;
% 
% missing_values.pd(np_test.np4.PD.idx_nan, 5) = 1;
% missing_values.hc(np_test.np4.HC.idx_nan, 5) = 1;
% 
% missing_values.pd(find(data_pd.GENETICS == "NA"), 6) = 1;
% missing_values.hc(find(data_hc.GENETICS == "NA"), 6) = 1;
% 
% missing_values.pd(find(data_pd.ANYFAMPD == "NA"), 7) = 1;
% missing_values.hc(find(data_hc.ANYFAMPD == "NA"), 7) = 1;
% 
% missing_values.pd(find(data_pd.ETHNICITY == "NA"), 8) = 1;
% missing_values.hc(find(data_hc.ETHNICITY == "NA"), 8) = 1;
% 
% missing_values.pd(find(data_pd.SEX == "NA"), 9) = 1;
% missing_values.hc(find(data_hc.SEX == "NA"), 9) = 1;
% 
% missing_values.pd(find(isnan(data_pd.ENROLL_AGE)), 10) = 1;
% missing_values.hc(find(isnan(data_hc.ENROLL_AGE)), 10) = 1;
% 
% missing_values.pd(find(isnan(data_pd.HTCM)), 11) = 1;
% missing_values.hc(find(isnan(data_hc.HTCM)), 11) = 1;
% 
% missing_values.pd(find(isnan(data_pd.WGTKG)), 12) = 1;
% missing_values.hc(find(isnan(data_hc.WGTKG)), 12) = 1;
% 
% 
% missing_values.pd(find(data_pd.HANDED == "NA"), 13) = 1;
% missing_values.hc(find(data_hc.HANDED == "NA"), 13) = 1;
% 
% missing_values.pd(find(string(data_pd.PRIMDIAG) == "NA"), 14) = 1;
% missing_values.hc(find(string(data_hc.PRIMDIAG) == "NA"), 14) = 1;
% 
% 
% figure
% subplot(121)
% imagesc(missing_values.hc)
% cmap = jet(2);
% colormap(cmap)
% hold on
% L = line(ones(2), ones(2), 'LineWidth', 2);
% set(L, {'color'}, mat2cell(cmap, ones(1, 2), 3));
% title('HC')
% xticks(1:14)
% xticklabels({'np1r', 'np1p', 'np2', 'np3', 'np4', 'genetics', 'familiarity', 'ethnicity', 'sex', 'age', 'height', 'weight', 'hand', 'primary diagnosis'})
% legend('Data', 'NaN')
% 
% subplot(122)
% imagesc(missing_values.pd)
% title('PD')
% cmap = jet(2);
% colormap(cmap)
% hold on
% L = line(ones(2), ones(2), 'LineWidth', 2);
% set(L, {'color'}, mat2cell(cmap, ones(1, 2), 3));
% xticks(1:14)
% xticklabels({'np1r', 'np1p', 'np2', 'np3', 'np4', 'genetics', 'familiarity', 'ethnicity', 'sex', 'age', 'height', 'weight', 'hand', 'primary diagnosis'})
% legend('Data', 'NaN')
% 
% 
% nan_count_pd = zeros(1, 14);
% nan_count_hc = zeros(1, 14);
% for i = 1:14
%     nan_count_pd(i) = sum(missing_values.pd(:,i));
%     nan_count_hc(i) = sum(missing_values.hc(:,i));
% end
% 
% figure
% subplot(121)
% bar(nan_count_hc, '')
% title('HC')
% xticks(1:14)
% xticklabels({'np1r', 'np1p', 'np2', 'np3', 'np4', 'genetics', 'familiarity', 'ethnicity', 'sex', 'age', 'height', 'weight', 'hand', 'primary diagnosis'})
% subplot(122)
% bar(nan_count_pd)
% title('PD')
% xticks(1:14)
% xticklabels({'np1r', 'np1p', 'np2', 'np3', 'np4', 'genetics', 'familiarity', 'ethnicity', 'sex', 'age', 'height', 'weight', 'hand', 'primary diagnosis'})
% 

%% ANALYSIS OF DEMOGRAPHICS DATA

%% CONTINUOUS VARIABLES

% age: already filled up the NaN with the date of the dat scan and the
% birth date
% height: PD: 23 NaN HC: 1 (find(isnan(variables.height.pd)))
% weight: PD: 21 NaN HC: 1(find(isnan(variables.weight.pd)))

% AGE --------------------------------------------------------------------
% Check gaussianity of AGE
[gaussianity.pd.age.h, gaussianity.pd.age.p] = lillietest(variables.age.pd);
if gaussianity.pd.age.h == 1
    disp('AGE in PD is NOT normally distributed')
else 
    disp('AGE in PD is normally distributed')
end

[gaussianity.hc.age.h, gaussianity.hc.age.p] = lillietest(variables.age.hc);
if gaussianity.hc.age.h == 1
    disp('AGE in HC is NOT normally distributed')
else
    disp('AGE in HC is normally distributed')
end

% Anova test AGE between PD and HC
y = [variables.age.pd; variables.age.hc]';
groups = [ones(1,size(variables.age.pd,1)), 2*ones(1,size(variables.age.hc,1))];
[anova_test.age.p,tbl,stats] = anova1(y, groups);
if anova_test.age.p > 0.05
    disp('AGE: Accepted null hyp (same mean)')
else 
    disp('AGE: Rejected null hyp (different mean)')
end

% HEIGHT -----------------------------------------------------------------
% Check gaussianity of HEIGHT
[gaussianity.pd.height.h, gaussianity.pd.height.p] = lillietest(variables.height.pd);
if gaussianity.pd.height.h == 1
    disp('HEIGHT in PD is NOT normally distributed')
else 
    disp('HEIGHT in PD is normally distributed')
end

[gaussianity.hc.height.h, gaussianity.hc.height.p] = lillietest(variables.height.hc);
if gaussianity.hc.height.h == 1
    disp('HEIGHT in HC is NOT normally distributed')
else
    disp('HEIGHT in HC is normally distributed')
end

% Anova test HEIGHT between PD and HC
y = [variables.height.pd; variables.height.hc]';
groups = [ones(1,size(variables.height.pd,1)), 2*ones(1,size(variables.height.hc,1))];
[anova_test.height.p,tbl,stats] = anova1(y, groups);
if anova_test.height.p > 0.05
    disp('HEIGHT: Accepted null hyp (same mean)')
else 
    disp('HEIGHT: Rejected null hyp (different mean)')
end

% WEIGHT -----------------------------------------------------------------
% Check gaussianity of WEIGHT
[gaussianity.pd.weight.h, gaussianity.pd.weight.p] = lillietest(variables.weight.pd);
if gaussianity.pd.weight.h == 1
    disp('WEIGHT in PD is NOT normally distributed')
else 
    disp('WEIGHT in PD is normally distributed')
end

[gaussianity.hc.weight.h, gaussianity.hc.weight.p] = lillietest(variables.weight.hc);
if gaussianity.hc.weight.h == 1
    disp('WEIGHT in HC is NOT normally distributed')
else
    disp('WEIGHT in HC is normally distributed')
end

% Anova test WEIGHT between PD and HC
y = [variables.weight.pd; variables.weight.hc]';
groups = [ones(1,size(variables.weight.pd,1)), 2*ones(1,size(variables.weight.hc,1))];
[anova_test.weight.p,tbl,stats] = anova1(y, groups);
if anova_test.weight.p > 0.05
    disp('WEIGHT: Accepted null hyp (same mean)')
else 
    disp('WEIGHT: Rejected null hyp (different mean)')
end

%% DISCRETE VARIABLES
% familiarity, ethnicity, sex, hand, genetics (solo PD)
% The Wilcoxon test function already removes the nan from the data

% FAMILIARITY -----------------------------------------------------------
% transform in numerical indeces
temp_pd = zeros(length(variables.familiarity.pd),1);
for i = 1:length(variables.familiarity.pd)
    if cell2mat(variables.familiarity.pd(i)) == 'NA'
        temp_pd(i) = NaN;
    elseif cell2mat(variables.familiarity.pd(i)) == '1'
        temp_pd(i) = 1;
    elseif cell2mat(variables.familiarity.pd(i)) == '0'
        temp_pd(i) = 0;
    end
end
temp_hc = zeros(length(variables.familiarity.hc),1);
for i = 1:length(variables.familiarity.hc)
    if cell2mat(variables.familiarity.hc(i)) == 'NA'
        temp_hc(i) = NaN;
    elseif cell2mat(variables.familiarity.hc(i)) == '1'
        temp_hc(i) = 1;
    elseif cell2mat(variables.familiarity.hc(i)) == '0'
        temp_hc(i) = 0;
    end
end
% save back in variables struct
variables.familiarity.pd = temp_pd;
variables.familiarity.hc = temp_hc;

% Perfom Wilcoxon test
[wilcoxon.familiarity.p,wilcoxon.familiarity.h,stats] = ranksum(variables.familiarity.pd,variables.familiarity.hc);
if wilcoxon.familiarity.h == 1
    disp('FAMILIARITY: Rejected null hyp (different mean)')
else
    disp('FAMILIARITY: Accepted null hyp (same mean)') 
end

clear temp_pd temp_hc

% ETHNICITY -----------------------------------------------------------
% transform in numerical indeces
values_ethnicity = unique(variables.ethnicity.pd);
temp_pd = zeros(length(variables.ethnicity.pd),1);
for i = 1:length(variables.ethnicity.pd)
    for j = 1:length(values_ethnicity)
        if string(cell2mat(variables.ethnicity.pd(i))) == string(cell2mat(values_ethnicity(j)))
            temp_pd(i) = j;
        end
    end
end
temp_hc = zeros(length(variables.ethnicity.hc),1);
for i = 1:length(variables.ethnicity.hc)
    for j = 1:length(values_ethnicity)
        if string(cell2mat(variables.ethnicity.hc(i))) == string(cell2mat(values_ethnicity(j)))
            temp_hc(i) = j;
        end
    end
end
% save back in variables struct
variables.ethnicity.pd = temp_pd;
variables.ethnicity.hc = temp_hc;

% Perfom Wilcoxon test
[wilcoxon.ethnicity.p,wilcoxon.ethnicity.h,stats] = ranksum(variables.ethnicity.pd,variables.ethnicity.hc);
if wilcoxon.ethnicity.h == 1
    disp('ETHNICITY: Rejected null hyp (different mean)')
else
    disp('ETHNICITY: Accepted null hyp (same mean)') 
end

clear temp_pd temp_hc

% SEX -----------------------------------------------------------
% transform in numerical indeces
values_sex = unique(variables.sex.pd);
temp_pd = zeros(length(variables.sex.pd),1);
for i = 1:length(variables.sex.pd)
    for j = 1:length(values_sex)
        if string(cell2mat(variables.sex.pd(i))) == string(cell2mat(values_sex(j)))
            temp_pd(i) = j;
        end
    end
end
temp_hc = zeros(length(variables.sex.hc),1);
for i = 1:length(variables.sex.hc)
    for j = 1:length(values_sex)
        if string(cell2mat(variables.sex.hc(i))) == string(cell2mat(values_sex(j)))
            temp_hc(i) = j;
        end
    end
end
% save back in variables struct
variables.sex.pd = temp_pd;
variables.sex.hc = temp_hc;

% Perfom Wilcoxon test
[wilcoxon.sex.p,wilcoxon.sex.h,stats] = ranksum(variables.sex.pd,variables.sex.hc);
if wilcoxon.sex.h == 1
    disp('SEX: Rejected null hyp (different mean)')
else
    disp('SEX: Accepted null hyp (same mean)') 
end

clear temp_pd temp_hc

% HAND -----------------------------------------------------------
% transform in numerical indeces
values_hand = unique(variables.hand.pd);
temp_pd = zeros(length(variables.hand.pd),1);
for i = 1:length(variables.hand.pd)
    for j = 1:length(values_hand)
        if string(cell2mat(variables.hand.pd(i))) == string(cell2mat(values_hand(j)))
            temp_pd(i) = j;
        end
        if string(cell2mat(values_hand(j))) == 'NA'
            temp_pd(i) = NaN;
        end
    end
end
temp_hc = zeros(length(variables.hand.hc),1);
for i = 1:length(variables.hand.hc)
    for j = 1:length(values_hand)
        if string(cell2mat(variables.hand.hc(i))) == string(cell2mat(values_hand(j)))
            temp_hc(i) = j;
        end
        if string(cell2mat(values_hand(j))) == 'NA'
            temp_hc(i) = NaN;
        end
    end
end
% save back in variables struct
variables.hand.pd = temp_pd;
variables.hand.hc = temp_hc;

% Perfom Wilcoxon test
[wilcoxon.hand.p,wilcoxon.hand.h,stats] = ranksum(variables.hand.pd,variables.hand.hc);
if wilcoxon.hand.h == 1
    disp('HAND: Rejected null hyp (different mean)')
else
    disp('HAND: Accepted null hyp (same mean)') 
end

clear temp_pd temp_hc







%% STATISTICAL ANALYSIS
%% ONE WAY ANOVA
% Control gauss
rois_names = fieldnames(DATSCAN);
for i=2:length(rois_names)
    rois_names_diag = fieldnames(DATSCAN.(rois_names{i,1}));
    for j =1:length(rois_names_diag)
        lil = lillietest(DATSCAN.(rois_names{i}).(rois_names_diag{j}));
        if lil == 1
            disp([rois_names{i}; rois_names_diag{j}; "not guassian"])
        end
    end
end
%%Absulte value
% Caudate
y = [abs(DATSCAN.CAUDATE_lat.PD)' abs(DATSCAN.CAUDATE_lat.HC)']';
group1 = [repmat("PD_caudate",length(DATSCAN.CAUDATE_lat.PD),1);repmat("HC_caudate",length(DATSCAN.CAUDATE_lat.HC),1)];
% [p,tbl,stats]  = anova1(y,group1);
% 
% % Putamen
y = [abs(DATSCAN.PUTAMEN_lat.PD)' abs(DATSCAN.PUTAMEN_lat.HC)']';
group2 = [repmat("PD_putamen",length(DATSCAN.PUTAMEN_lat.PD),1);repmat("HC_putamen",length(DATSCAN.PUTAMEN_lat.HC),1)];
[p,tbl,stats]  = anova1(y,group2);
% 
% % Putamen ANT
y = [abs(DATSCAN.PUTAMEN_ANT_lat.PD)' abs(DATSCAN.PUTAMEN_ANT_lat.HC)']';
group3 = [repmat("PD_putamen_ant",length(DATSCAN.PUTAMEN_ANT_lat.PD),1);repmat("HC_putamen_ant",length(DATSCAN.PUTAMEN_ANT_lat.HC),1)];
[p,tbl,stats]  = anova1(y,group3);
% 
y = [abs(DATSCAN.CAUDATE_lat.PD)' abs(DATSCAN.CAUDATE_lat.HC)' abs(DATSCAN.PUTAMEN_lat.PD)' abs(DATSCAN.PUTAMEN_lat.HC)' abs(DATSCAN.PUTAMEN_ANT_lat.PD)'  abs(DATSCAN.PUTAMEN_ANT_lat.HC)']';
group4 = [group1; group2; group3];
[p,tbl,stats]  = anova1(y,group4);
% 

%multcompare(stats)
%%Left
%abs(DATSCAN.PUTAMEN_lat.PD)' abs(DATSCAN.PUTAMEN_lat.HC)' abs(DATSCAN.PUTAMEN_ANT_lat.PD)'  abs(DATSCAN.PUTAMEN_ANT_lat.HC)'
% y = [DATSCAN.CAUDATE_lat.PD(lateralization.CAUDATE.PD.left.index)' DATSCAN.CAUDATE_lat.HC(lateralization.CAUDATE.HC.left.index)']';
% group5 = [repmat("PD_caudate_LEFT",length(lateralization.CAUDATE.PD.left.index),1);repmat("HC_caudate_LEFT",length(lateralization.CAUDATE.HC.left.index),1)];
% [p,tbl,stats]  = anova1(y,group5);

%% CORRELATION MATRIX

%% all important variables ---------------------------------------------
new_data_hc(:, "CAUDATE LAT") = table(DATSCAN.CAUDATE_lat.HC(:));
new_data_hc(:, "PUTAMEN LAT") = table(DATSCAN.PUTAMEN_lat.HC(:));
new_data_hc(:, "PUTAMEN ANT LAT") = table(DATSCAN.PUTAMEN_ANT_lat.HC(:));

new_data_pd(:, "CAUDATE LAT") = table(DATSCAN.CAUDATE_lat.PD(:));
new_data_pd(:, "PUTAMEN LAT") = table(DATSCAN.PUTAMEN_lat.PD(:));
new_data_pd(:, "PUTAMEN ANT LAT") = table(DATSCAN.PUTAMEN_ANT_lat.PD(:));

col_to_keep = [27:55, 61:94];

new_data_hc = [new_data_hc data_hc(:, col_to_keep)];
new_data_pd = [new_data_pd data_pd(:, col_to_keep)];

% CORRELATION MATRIX
R_HC = corrcoef(table2array(new_data_hc), 'Rows', 'complete');
R_PD = corrcoef(table2array(new_data_pd), 'Rows', 'complete');

figure
imagesc(R_HC)
colormap parula 
colorbar
xticks(1:70)
xticklabels(new_data_hc.Properties.VariableNames)
yticks(1:70)
yticklabels(new_data_hc.Properties.VariableNames)
title("Correlation for HC")

figure
imagesc(R_PD)
colormap parula
colorbar
xticks(1:70)
xticklabels(new_data_hc.Properties.VariableNames)
yticks(1:70)
yticklabels(new_data_hc.Properties.VariableNames)
title("Correlation for PD")

%% correlation between lateralization and np3 test -----------------------
new_data_hc_np3 = new_data_hc(:, [1:3, 33:end]);
new_data_pd_np3 = new_data_pd(:, [1:3, 33:end]);

% CORRELATION MATRIX
R_np3_HC = corrcoef(table2array(new_data_hc_np3), 'Rows', 'complete');
R_np3_PD = corrcoef(table2array(new_data_pd_np3), 'Rows', 'complete');

figure
imagesc(R_np3_HC)
colormap parula 
colorbar
xticks(1:width(new_data_hc_np3))
xticklabels(new_data_hc_np3.Properties.VariableNames)
yticks(1:width(new_data_hc_np3))
yticklabels(new_data_hc_np3.Properties.VariableNames)
title("Correlation for HC - NP3 TEST")

figure
imagesc(R_np3_PD)
colormap parula
colorbar
xticks(1:width(new_data_hc_np3))
xticklabels(new_data_hc_np3.Properties.VariableNames)
yticks(1:width(new_data_hc_np3))
yticklabels(new_data_hc_np3.Properties.VariableNames)
title("Correlation for PD - NP3 TEST")

%% correlation between CAUDATE RIGHT lateralization and np3 test RIGHT
% non ha senso per HC che non ha lateralizzazione, quindi solo PD

np3_right = [1, 7, 9, 11, 13, 15, 17, 19, 27, 29, 31, 33];
caudate_new_data_pd_np3_rr = new_data_pd_np3(lateralization.CAUDATE.PD.right.index,np3_right);

% CORRELATION MATRIX
R_caudate_np3_PD_rr = corrcoef(table2array(caudate_new_data_pd_np3_rr), 'Rows', 'complete');

figure
imagesc(R_caudate_np3_PD_rr)
colormap parula
colorbar
xticks(1:width(caudate_new_data_pd_np3_rr))
xticklabels(caudate_new_data_pd_np3_rr.Properties.VariableNames)
yticks(1:width(caudate_new_data_pd_np3_rr))
yticklabels(caudate_new_data_pd_np3_rr.Properties.VariableNames)
title("Correlation PD CAUDATE lat RIGHT - NP3 TEST RIGHT")

% %% correlation between CAUDATE RIGHT lateralization and np3 test LEFT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_left = [1, 8, 10, 12, 14, 16, 18, 20, 28, 30, 32, 34];
% caudate_new_data_pd_np3_rl = new_data_pd_np3(lateralization.CAUDATE.PD.right.index,np3_left);
% 
% % CORRELATION MATRIX
% R_caudate_np3_PD_rl = corrcoef(table2array(caudate_new_data_pd_np3_rl), 'Rows', 'complete');
% 
% figure
% imagesc(R_caudate_np3_PD_rl)
% colormap parula
% colorbar
% xticks(1:width(caudate_new_data_pd_np3_rl))
% xticklabels(caudate_new_data_pd_np3_rl.Properties.VariableNames)
% yticks(1:width(caudate_new_data_pd_np3_rl))
% yticklabels(caudate_new_data_pd_np3_rl.Properties.VariableNames)
% title("Correlation PD CAUDATE lat RIGHT - NP3 TEST LEFT")
% 
% %% correlation between CAUDATE LEFT lateralization and np3 test RIGHT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_right = [1, 7, 9, 11, 13, 15, 17, 19, 27, 29, 31, 33];
% caudate_new_data_pd_np3_lr = new_data_pd_np3(lateralization.CAUDATE.PD.left.index,np3_right);
% 
% % CORRELATION MATRIX
% R_caudate_np3_PD_lr = corrcoef(table2array(caudate_new_data_pd_np3_lr), 'Rows', 'complete');
% 
% figure
% imagesc(R_caudate_np3_PD_lr)
% colormap parula
% colorbar
% xticks(1:width(caudate_new_data_pd_np3_lr))
% xticklabels(caudate_new_data_pd_np3_lr.Properties.VariableNames)
% yticks(1:width(caudate_new_data_pd_np3_lr))
% yticklabels(caudate_new_data_pd_np3_lr.Properties.VariableNames)
% title("Correlation PD CAUDATE lat LEFT - NP3 TEST RIGHT")
% 
% %% correlation between CAUDATE LEFT lateralization and np3 test LEFT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_left = [1, 8, 10, 12, 14, 16, 18, 20, 28, 30, 32, 34];
% caudate_new_data_pd_np3_ll = new_data_pd_np3(lateralization.CAUDATE.PD.left.index,np3_left);
% 
% % CORRELATION MATRIX
% R_caudate_np3_PD_ll = corrcoef(table2array(caudate_new_data_pd_np3_ll), 'Rows', 'complete');
% 
% figure
% imagesc(R_caudate_np3_PD_ll)
% colormap parula
% colorbar
% xticks(1:width(caudate_new_data_pd_np3_ll))
% xticklabels(caudate_new_data_pd_np3_ll.Properties.VariableNames)
% yticks(1:width(caudate_new_data_pd_np3_ll))
% yticklabels(caudate_new_data_pd_np3_ll.Properties.VariableNames)
% title("Correlation PD CAUDATE lat LEFT - NP3 TEST LEFT")
% 
% %% correlation between PUTAMEN RIGHT lateralization and np3 test RIGHT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_right = [2, 7, 9, 11, 13, 15, 17, 19, 27, 29, 31, 33];
% putamen_new_data_pd_np3_rr = new_data_pd_np3(lateralization.PUTAMEN.PD.right.index,np3_right);
% 
% % CORRELATION MATRIX
% R_putamen_np3_PD_rr = corrcoef(table2array(putamen_new_data_pd_np3_rr), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_np3_PD_rr)
% colormap parula
% colorbar
% xticks(1:width(putamen_new_data_pd_np3_rr))
% xticklabels(putamen_new_data_pd_np3_rr.Properties.VariableNames)
% yticks(1:width(putamen_new_data_pd_np3_rr))
% yticklabels(putamen_new_data_pd_np3_rr.Properties.VariableNames)
% title("Correlation PD PUTAMEN lat RIGHT - NP3 TEST RIGHT")
% 
% %% correlation between PUTAMEN RIGHT lateralization and np3 test LEFT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_left = [2, 8, 10, 12, 14, 16, 18, 20, 28, 30, 32, 34];
% putamen_new_data_pd_np3_rl = new_data_pd_np3(lateralization.PUTAMEN.PD.right.index,np3_left);
% 
% % CORRELATION MATRIX
% R_putamen_np3_PD_rl = corrcoef(table2array(putamen_new_data_pd_np3_rl), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_np3_PD_rl)
% colormap parula
% colorbar
% xticks(1:width(putamen_new_data_pd_np3_rl))
% xticklabels(putamen_new_data_pd_np3_rl.Properties.VariableNames)
% yticks(1:width(putamen_new_data_pd_np3_rl))
% yticklabels(putamen_new_data_pd_np3_rl.Properties.VariableNames)
% title("Correlation PD PUTAMEN lat RIGHT - NP3 TEST LEFT")
% 
% %% correlation between PUTAMEN LEFT lateralization and np3 test RIGHT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_right = [2, 7, 9, 11, 13, 15, 17, 19, 27, 29, 31, 33];
% putamen_new_data_pd_np3_lr = new_data_pd_np3(lateralization.PUTAMEN.PD.left.index,np3_right);
% 
% % CORRELATION MATRIX
% R_putamen_np3_PD_lr = corrcoef(table2array(putamen_new_data_pd_np3_lr), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_np3_PD_lr)
% colormap parula
% colorbar
% xticks(1:width(putamen_new_data_pd_np3_lr))
% xticklabels(putamen_new_data_pd_np3_lr.Properties.VariableNames)
% yticks(1:width(putamen_new_data_pd_np3_lr))
% yticklabels(putamen_new_data_pd_np3_lr.Properties.VariableNames)
% title("Correlation PD PUTAMEN lat LEFT - NP3 TEST RIGHT")
% 
% %% correlation between PUTAMEN LEFT lateralization and np3 test LEFT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_left = [2, 8, 10, 12, 14, 16, 18, 20, 28, 30, 32, 34];
% putamen_new_data_pd_np3_ll = new_data_pd_np3(lateralization.PUTAMEN.PD.left.index,np3_left);
% 
% % CORRELATION MATRIX
% R_putamen_np3_PD_ll = corrcoef(table2array(putamen_new_data_pd_np3_ll), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_np3_PD_ll)
% colormap parula
% colorbar
% xticks(1:width(putamen_new_data_pd_np3_ll))
% xticklabels(putamen_new_data_pd_np3_ll.Properties.VariableNames)
% yticks(1:width(putamen_new_data_pd_np3_ll))
% yticklabels(putamen_new_data_pd_np3_ll.Properties.VariableNames)
% title("Correlation PD PUTAMEN lat LEFT - NP3 TEST LEFT")
% 
% %% correlation between PUTAMEN ANT RIGHT lateralization and np3 test RIGHT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_right = [3, 7, 9, 11, 13, 15, 17, 19, 27, 29, 31, 33];
% putamen_ant_new_data_pd_np3_rr = new_data_pd_np3(lateralization.PUTAMEN_ANT.PD.right.index,np3_right);
% 
% % CORRELATION MATRIX
% R_putamen_ant_np3_PD_rr = corrcoef(table2array(putamen_ant_new_data_pd_np3_rr), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_ant_np3_PD_rr)
% colormap parula
% colorbar
% xticks(1:width(putamen_ant_new_data_pd_np3_rr))
% xticklabels(putamen_ant_new_data_pd_np3_rr.Properties.VariableNames)
% yticks(1:width(putamen_ant_new_data_pd_np3_rr))
% yticklabels(putamen_ant_new_data_pd_np3_rr.Properties.VariableNames)
% title("Correlation PD PUTAMEN ANT lat RIGHT - NP3 TEST RIGHT")
% 
% %% correlation between PUTAMEN RIGHT lateralization and np3 test LEFT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_left = [3, 8, 10, 12, 14, 16, 18, 20, 28, 30, 32, 34];
% putamen_ant_new_data_pd_np3_rl = new_data_pd_np3(lateralization.PUTAMEN_ANT.PD.right.index,np3_left);
% 
% % CORRELATION MATRIX
% R_putamen_ant_np3_PD_rl = corrcoef(table2array(putamen_ant_new_data_pd_np3_rl), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_ant_np3_PD_rl)
% colormap parula
% colorbar
% xticks(1:width(putamen_ant_new_data_pd_np3_rl))
% xticklabels(putamen_ant_new_data_pd_np3_rl.Properties.VariableNames)
% yticks(1:width(putamen_ant_new_data_pd_np3_rl))
% yticklabels(putamen_ant_new_data_pd_np3_rl.Properties.VariableNames)
% title("Correlation PD PUTAMEN ANT lat RIGHT - NP3 TEST LEFT")
% 
% %% correlation between PUTAMEN LEFT lateralization and np3 test RIGHT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_right = [3, 7, 9, 11, 13, 15, 17, 19, 27, 29, 31, 33];
% putamen_ant_new_data_pd_np3_lr = new_data_pd_np3(lateralization.PUTAMEN_ANT.PD.left.index,np3_right);
% 
% % CORRELATION MATRIX
% R_putamen_ant_np3_PD_lr = corrcoef(table2array(putamen_ant_new_data_pd_np3_lr), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_ant_np3_PD_lr)
% colormap parula
% colorbar
% xticks(1:width(putamen_ant_new_data_pd_np3_lr))
% xticklabels(putamen_ant_new_data_pd_np3_lr.Properties.VariableNames)
% yticks(1:width(putamen_ant_new_data_pd_np3_lr))
% yticklabels(putamen_ant_new_data_pd_np3_lr.Properties.VariableNames)
% title("Correlation PD PUTAMEN ANT lat LEFT - NP3 TEST RIGHT")
% 
% %% correlation between PUTAMEN ANT LEFT lateralization and np3 test LEFT
% % non ha senso per HC che non ha lateralizzazione, quindi solo PD
% 
% np3_left = [3, 8, 10, 12, 14, 16, 18, 20, 28, 30, 32, 34];
% putamen_ant_new_data_pd_np3_ll = new_data_pd_np3(lateralization.PUTAMEN_ANT.PD.left.index,np3_left);
% 
% % CORRELATION MATRIX
% R_putamen_ant_np3_PD_ll = corrcoef(table2array(putamen_ant_new_data_pd_np3_ll), 'Rows', 'complete');
% 
% figure
% imagesc(R_putamen_ant_np3_PD_ll)
% colormap parula
% colorbar
% xticks(1:width(putamen_ant_new_data_pd_np3_ll))
% xticklabels(putamen_ant_new_data_pd_np3_ll.Properties.VariableNames)
% yticks(1:width(putamen_ant_new_data_pd_np3_ll))
% yticklabels(putamen_ant_new_data_pd_np3_ll.Properties.VariableNames)
% title("Correlation PD PUTAMEN ANT lat LEFT - NP3 TEST LEFT")

%% LATERALIZATION AND DOMINANT HAND
% count_pd_latR_handR = length()
%data_pd_lat_dx_mano_dx



%% SCATTERPLOT
j = [11,18,32,66];
for j=j
    figure
    for i=1:3
        subplot(1,3,i)
        scatter(table2array(new_data_pd(:,i)),table2array(new_data_pd(:,j)))
        hold on
        scatter(table2array(new_data_hc(:,i)), table2array(new_data_hc(:,j)))
        xlabel('Lateralization')
        ylabel('Symptoms')
        legend('PD','HC')
        hold off
        title(new_data_pd.Properties.VariableNames{j}, new_data_pd.Properties.VariableNames{i})
    end
end

%% LINEAR REGRESSION of variables of interest
%% right
figure
subplot(131)
symptoms_r = [data_pd.NP3RIGN,data_pd.NP3RIGRU,data_pd.NP3RIGRL,data_pd.NP3PTRMR,data_pd.NP3KTRMR];
model_caud_r = fitlm(symptoms_r,DATSCAN.CAUDATE_lat.PD);
plot(model_caud_r)
subplot(132)

symptoms_r = [data_pd.NP3RIGN,data_pd.NP3RIGRU,data_pd.NP3RIGRL,data_pd.NP3PTRMR,data_pd.NP3KTRMR];
model_put_ant_r = fitlm(symptoms_r,DATSCAN.PUTAMEN_ANT_lat.PD);
plot(model_put_ant_r)
subplot(133)

symptoms_r = [data_pd.NP3RIGN,data_pd.NP3RIGRU,data_pd.NP3RIGRL,data_pd.NP3PTRMR,data_pd.NP3KTRMR];
model_put_r = fitlm(symptoms_r,DATSCAN.PUTAMEN_lat.PD);
plot(model_put_r)

% Statistics
%caudate
% anova
anova_caud = anova(model_caud_r);
% coef test
coed_test_caud = coefTest(model_caud_r);
% [pdep_caud,x,y] = partialDependence(model_caud,{'x1','x3'});
% figure
% imagesc(x,y,pdep_caud)
% putamen
% anova
anova_put = anova(model_put_r);
% coef test
% coed_test_caud = coefTest(model_put);
% [pdep_put,x,y] = partialDependence(model_put,{'x1','x3'});
% figure
% imagesc(x,y,pdep_put)
%% left
figure
subplot(131)
symptoms_l = [data_pd.NP3RIGN,data_pd.NP3RIGLU,data_pd.NP3RIGLL,data_pd.NP3PTRML,data_pd.NP3KTRML];
model_caud_l = fitlm(symptoms_l,DATSCAN.CAUDATE_lat.PD);
plot(model_caud_l)

subplot(132)
model_put_ant_l = fitlm(symptoms_l,DATSCAN.PUTAMEN_ANT_lat.PD);
plot(model_put_ant_l)

subplot(133)
model_put_l = fitlm(symptoms_l,DATSCAN.PUTAMEN_lat.PD);
plot(model_put_l)

% Statistics
%caudate
% anova
anova_caud = anova(model_caud_l);
% coef test
coed_test_caud = coefTest(model_caud_l);
% [pdep_caud,x,y] = partialDependence(model_caud,{'x1','x3'});
% figure
% imagesc(x,y,pdep_caud)
% putamen
% anova
anova_put = anova(model_put_l);
% coef test
% coed_test_caud = coefTest(model_put);
% [pdep_put,x,y] = partialDependence(model_put,{'x1','x3'});
% figure
% imagesc(x,y,pdep_put)
%% Control if lateralization in one direction is constant in all ROIs
%%PD
pd_patient_id_right_lat_caudate = lateralization.CAUDATE.PD.right.population.Var1;
pd_patient_id_right_lat_putamen = lateralization.PUTAMEN.PD.right.population.Var1;
pd_patient_id_right_lat_putamen_ant = lateralization.PUTAMEN_ANT.PD.right.population.Var1;

pd_patient_id_left_lat_caudate = lateralization.CAUDATE.PD.left.population.Var1;
pd_patient_id_left_lat_putamen = lateralization.PUTAMEN.PD.left.population.Var1;
pd_patient_id_left_lat_putamen_ant = lateralization.PUTAMEN_ANT.PD.left.population.Var1;

pd_patient_id_none_lat_caudate = lateralization.CAUDATE.PD.none.population.Var1;
pd_patient_id_none_lat_putamen = lateralization.PUTAMEN.PD.none.population.Var1;
pd_patient_id_none_lat_putamen_ant = lateralization.PUTAMEN_ANT.PD.none.population.Var1;


% Right lat
pd_caudate_putamen_right = intersect(pd_patient_id_right_lat_caudate,pd_patient_id_right_lat_putamen);
pd_caudate_putamen_ant_right = intersect(pd_patient_id_right_lat_caudate,pd_patient_id_right_lat_putamen_ant);
pd_putamen_putamen_ant_right = intersect(pd_patient_id_right_lat_putamen,pd_patient_id_right_lat_putamen_ant);
pd_caudate_putamen_putamen_ant_right = intersect(pd_caudate_putamen_right,pd_caudate_putamen_ant_right);

pd_caudate_only_right = setdiff(setdiff(pd_patient_id_right_lat_caudate,pd_caudate_putamen_ant_right),pd_caudate_putamen_right);
pd_putamen_only_right = setdiff(setdiff(pd_patient_id_right_lat_putamen,pd_caudate_putamen_ant_right),pd_putamen_putamen_ant_right);
pd_putamen_ant_only_right = setdiff(setdiff(pd_patient_id_right_lat_putamen_ant,pd_caudate_putamen_ant_right),pd_putamen_putamen_ant_right);

pd_caudate_right_putamen_ant_none = intersect(pd_patient_id_right_lat_caudate,pd_patient_id_none_lat_putamen_ant);
pd_caudate_right_putamen_none = intersect(pd_patient_id_right_lat_caudate,pd_patient_id_none_lat_putamen);
pd_caudate_right_others_none = intersect(pd_caudate_right_putamen_ant_none,pd_caudate_right_putamen_none);

pd_putamen_right_putamen_ant_none = intersect(pd_patient_id_right_lat_putamen,pd_patient_id_none_lat_putamen_ant);
pd_putamen_right_caudate_none = intersect(pd_patient_id_right_lat_putamen,pd_patient_id_none_lat_caudate);
pd_putamen_right_others_none = intersect(pd_putamen_right_caudate_none,pd_putamen_right_putamen_ant_none);

pd_putamen_ant_right_putamen_none = intersect(pd_patient_id_right_lat_putamen_ant,pd_patient_id_none_lat_putamen);
pd_putamen_ant_right_caudate_none = intersect(pd_patient_id_right_lat_putamen_ant,pd_patient_id_none_lat_caudate);
pd_putamen_ant_right_others_none = intersect(pd_putamen_ant_right_putamen_none,pd_putamen_ant_right_caudate_none);


% Left
pd_caudate_putamen_left = intersect(pd_patient_id_left_lat_caudate,pd_patient_id_left_lat_putamen);
pd_caudate_putamen_ant_left = intersect(pd_patient_id_left_lat_caudate,pd_patient_id_left_lat_putamen_ant);
pd_putamen_putamen_ant_left = intersect(pd_patient_id_left_lat_putamen,pd_patient_id_left_lat_putamen_ant);
pd_caudate_putamen_putamen_ant_left = intersect(pd_caudate_putamen_left,pd_caudate_putamen_ant_left);

pd_caudate_only_left = setdiff(setdiff(pd_patient_id_left_lat_caudate,pd_caudate_putamen_ant_left),pd_caudate_putamen_left);
pd_putamen_only_left = setdiff(setdiff(pd_patient_id_left_lat_putamen,pd_caudate_putamen_ant_left),pd_putamen_putamen_ant_left);
pd_putamen_ant_only_left = setdiff(setdiff(pd_patient_id_left_lat_putamen_ant,pd_caudate_putamen_ant_left),pd_putamen_putamen_ant_left);

pd_caudate_left_putamen_none = intersect(pd_patient_id_left_lat_caudate,pd_patient_id_none_lat_putamen);
pd_caudate_left_putamen_ant_none = intersect(pd_patient_id_left_lat_caudate,pd_patient_id_none_lat_putamen_ant);
pd_caudate_left_others_none = intersect(pd_caudate_left_putamen_none,pd_caudate_left_putamen_ant_none);

pd_putamen_left_putamen_ant_none = intersect(pd_patient_id_left_lat_putamen,pd_patient_id_none_lat_putamen_ant);
pd_putamen_left_caudate_none = intersect(pd_patient_id_left_lat_putamen,pd_patient_id_none_lat_caudate);
pd_putamen_left_others_none = intersect(pd_putamen_left_putamen_ant_none,pd_putamen_left_caudate_none);

pd_putamen_ant_left_putamen_none = intersect(pd_patient_id_left_lat_putamen_ant,pd_patient_id_none_lat_putamen);
pd_putamen_ant_left_caudate_none = intersect(pd_patient_id_left_lat_putamen_ant,pd_patient_id_none_lat_caudate);
pd_putamen_ant_left_others_none = intersect(pd_putamen_ant_left_putamen_none,pd_putamen_ant_left_caudate_none);


% Right and left
pd_caudate_right_putamen_ant_left = intersect(pd_patient_id_right_lat_caudate,pd_patient_id_left_lat_putamen_ant);
pd_caudate_right_putamen_left = intersect(pd_patient_id_right_lat_caudate,pd_patient_id_left_lat_putamen);
pd_caudate_right_others_left = intersect(pd_caudate_right_putamen_ant_left,pd_caudate_right_putamen_left);

pd_putamen_right_putamen_ant_left = intersect(pd_patient_id_right_lat_putamen,pd_patient_id_left_lat_putamen_ant);
pd_putamen_right_caudate_left = intersect(pd_patient_id_right_lat_putamen,pd_patient_id_left_lat_caudate);
pd_putamen_right_others_left = intersect(pd_putamen_right_putamen_ant_left,pd_putamen_right_caudate_left);

pd_putamen_ant_right_putamen_left = intersect(pd_patient_id_right_lat_putamen_ant,pd_patient_id_left_lat_putamen);
pd_putamen_ant_right_caudate_left = intersect(pd_patient_id_right_lat_putamen_ant,pd_patient_id_left_lat_caudate);
pd_putamen_ant_right_others_left = intersect(pd_putamen_ant_right_putamen_left,pd_putamen_ant_right_caudate_left);

% %%HC
% hc_patient_id_right_lat_caudate = lateralization.CAUDATE.HC.right.population.Var1;
% hc_patient_id_right_lat_putamen = lateralization.PUTAMEN.HC.right.population.Var1;
% hc_patient_id_right_lat_putamen_ant = lateralization.PUTAMEN_ANT.HC.right.population.Var1;
% 
% hc_patient_id_left_lat_caudate = lateralization.CAUDATE.HC.left.population.Var1;
% hc_patient_id_left_lat_putamen = lateralization.PUTAMEN.HC.left.population.Var1;
% hc_patient_id_left_lat_putamen_ant = lateralization.PUTAMEN_ANT.HC.left.population.Var1;
% 
% hc_patient_id_none_lat_caudate = lateralization.CAUDATE.HC.none.population.Var1;
% hc_patient_id_none_lat_putamen = lateralization.PUTAMEN.HC.none.population.Var1;
% hc_patient_id_none_lat_putamen_ant = lateralization.PUTAMEN_ANT.HC.none.population.Var1;
% 
% 
% % Right 
% hc_caudate_putamen_right = intersect(hc_patient_id_right_lat_caudate,hc_patient_id_right_lat_putamen);
% hc_caudate_putamen_ant_right = intersect(hc_patient_id_right_lat_caudate,hc_patient_id_right_lat_putamen_ant);
% hc_putamen_putamen_ant_right = intersect(hc_patient_id_right_lat_putamen,hc_patient_id_right_lat_putamen_ant);
% hc_caudate_putamen_putamen_ant_right = intersect(hc_caudate_putamen_right,hc_caudate_putamen_ant_right);
% 
% hc_caudate_only_right = setdiff(setdiff(hc_patient_id_right_lat_caudate,hc_caudate_putamen_ant_right),hc_caudate_putamen_right);
% hc_putamen_only_right = setdiff(setdiff(hc_patient_id_right_lat_putamen,hc_caudate_putamen_ant_right),hc_putamen_putamen_ant_right);
% hc_putamen_ant_only_right = setdiff(setdiff(hc_patient_id_right_lat_putamen_ant,hc_caudate_putamen_ant_right),hc_putamen_putamen_ant_right);
% 
% hc_caudate_right_putamen_ant_none = intersect(hc_patient_id_right_lat_caudate,hc_patient_id_none_lat_putamen_ant);
% hc_caudate_right_putamen_none = intersect(hc_patient_id_right_lat_caudate,hc_patient_id_none_lat_putamen);
% hc_caudate_right_others_none = intersect(hc_caudate_right_putamen_ant_none,hc_caudate_right_putamen_none);
% 
% hc_putamen_right_putamen_ant_none = intersect(hc_patient_id_right_lat_putamen,hc_patient_id_none_lat_putamen_ant);
% hc_putamen_right_caudate_none = intersect(hc_patient_id_right_lat_putamen,hc_patient_id_none_lat_caudate);
% hc_putamen_right_others_none = intersect(hc_putamen_right_caudate_none,hc_putamen_right_putamen_ant_none);
% 
% hc_putamen_ant_right_putamen_none = intersect(hc_patient_id_right_lat_putamen_ant,hc_patient_id_none_lat_putamen);
% hc_putamen_ant_right_caudate_none = intersect(hc_patient_id_right_lat_putamen_ant,hc_patient_id_none_lat_caudate);
% hc_putamen_ant_right_others_none = intersect(hc_putamen_ant_right_putamen_none,hc_putamen_ant_right_caudate_none);
% 
% 
% % Left
% hc_caudate_putamen_left = intersect(hc_patient_id_left_lat_caudate,hc_patient_id_left_lat_putamen);
% hc_caudate_putamen_ant_left = intersect(hc_patient_id_left_lat_caudate,hc_patient_id_left_lat_putamen_ant);
% hc_putamen_putamen_ant_left = intersect(hc_patient_id_left_lat_putamen,hc_patient_id_left_lat_putamen_ant);
% hc_caudate_putamen_putamen_ant_left = intersect(hc_caudate_putamen_left,hc_caudate_putamen_ant_left);
% 
% hc_caudate_only_left = setdiff(setdiff(hc_patient_id_left_lat_caudate,hc_caudate_putamen_ant_left),hc_caudate_putamen_left);
% hc_putamen_only_left = setdiff(setdiff(hc_patient_id_left_lat_putamen,hc_caudate_putamen_ant_left),hc_putamen_putamen_ant_left);
% hc_putamen_ant_only_left = setdiff(setdiff(hc_patient_id_left_lat_putamen_ant,hc_caudate_putamen_ant_left),hc_putamen_putamen_ant_left);
% 
% hc_caudate_left_putamen_none = intersect(hc_patient_id_left_lat_caudate,hc_patient_id_none_lat_putamen);
% hc_caudate_left_putamen_ant_none = intersect(pd_patient_id_left_lat_caudate,hc_patient_id_none_lat_putamen_ant);
% hc_caudate_left_others_none = intersect(hc_caudate_left_putamen_none,hc_caudate_left_putamen_ant_none);
% 
% hc_putamen_left_putamen_ant_none = intersect(hc_patient_id_left_lat_putamen,hc_patient_id_none_lat_putamen_ant);
% hc_putamen_left_caudate_none = intersect(hc_patient_id_left_lat_putamen,hc_patient_id_none_lat_caudate);
% hc_putamen_left_others_none = intersect(hc_putamen_left_putamen_ant_none,hc_putamen_left_caudate_none);
% 
% hc_putamen_ant_left_putamen_none = intersect(hc_patient_id_left_lat_putamen_ant,hc_patient_id_none_lat_putamen);
% hc_putamen_ant_left_caudate_none = intersect(hc_patient_id_left_lat_putamen_ant,hc_patient_id_none_lat_caudate);
% hc_putamen_ant_left_others_none = intersect(hc_putamen_ant_left_putamen_none,hc_putamen_ant_left_caudate_none);
% 
% 
% % Right and left
% hc_caudate_right_putamen_ant_left = intersect(hc_patient_id_right_lat_caudate,hc_patient_id_left_lat_putamen_ant);
% hc_caudate_right_putamen_left = intersect(hc_patient_id_right_lat_caudate,hc_patient_id_left_lat_putamen);
% hc_caudate_right_others_left = intersect(hc_caudate_right_putamen_ant_left,hc_caudate_right_putamen_left);
% 
% hc_putamen_right_putamen_ant_left = intersect(hc_patient_id_right_lat_putamen,hc_patient_id_left_lat_putamen_ant);
% hc_putamen_right_caudate_left = intersect(hc_patient_id_right_lat_putamen,pd_patient_id_left_lat_caudate);
% hc_putamen_right_others_left = intersect(hc_putamen_right_putamen_ant_left,hc_putamen_right_caudate_left);
% 
% hc_putamen_ant_right_putamen_left = intersect(hc_patient_id_right_lat_putamen_ant,hc_patient_id_left_lat_putamen);
% hc_putamen_ant_right_caudate_left = intersect(hc_patient_id_right_lat_putamen_ant,pd_patient_id_left_lat_caudate);
% hc_putamen_ant_right_others_left = intersect(hc_putamen_ant_right_putamen_left,hc_putamen_ant_right_caudate_left);

% per HC solo due outliers con lat sx in putamen
%% LASSO FEATURE SELECTION



