%%% BIOMARKER, PRECISION MEDICINE AND DRUG DEVELOPMENT - HOMEWORK
%%% Benedetta Corso, Dalila Dattoli, Letizia Rossato - Group 12

% Investigation brain dopamine lateralization in Parkinson's Disease
clear all
close all
clc

%% load data
data = readtable("Patient_Master.csv");

% Regions of interest:
ROIs_labels = ["CAUDATE", "PUTAMEN", "PUTAMEN_ANT"];

% Check repetitions in PATNO - patients ID
IDs = data.PATNO;
if length(IDs) == length(unique(IDs))
    disp('No repetitions of patients IDs')
end
clear IDs

% Dataset sub-divisions
idx_samples.COHORT.HC = find(string(data.COHORT)=='HC'); % Healthy Controls
idx_samples.COHORT.PD = find(string(data.COHORT)=='PD'); % Parkinson's Disease
idx_samples.COHORT.SWEDD = find(string(data.COHORT)=='SWEDD'); % Scans without evidence of dopaminergic deficit
idx_samples.COHORT.Prodromal = find(string(data.COHORT)=='Prodromal'); % Early signs of PD - NO ONE

% work datasets
data_pd = data([idx_samples.COHORT.Prodromal; idx_samples.COHORT.PD; idx_samples.COHORT.SWEDD],:);
data_hc = data(idx_samples.COHORT.HC, :);

% SWEDD patients exclusion
ind_97_hc = find(data_hc.PRIMDIAG == 97);
ind_97_pd = find(data_pd.PRIMDIAG == 97);
for i = 1:length(ind_97_hc)
    data_hc(ind_97_hc(i),"PRIMDIAG") = {"NA"};
end
for i = 1:length(ind_97_pd)
    data_pd(ind_97_pd(i),"PRIMDIAG") = {"NA"};
end
clear ind_97_pd ind_97_hc


%% HANDLING MISSING VALUES

%%   Conversion to double after visual inspection
ind_to_convert = [11, 15, 23:24, 58, 96, 98:104, 118];
variables = data_pd.Properties.VariableNames;
for i = 1:length(ind_to_convert)
    data_pd.(variables{ind_to_convert(i)}) = str2double(data_pd.(variables{ind_to_convert(i)}));
    data_hc.(variables{ind_to_convert(i)}) = str2double(data_hc.(variables{ind_to_convert(i)}));
end
clear i ind_to_convert

%%   Extraction of NaNs indexes
% PD --------------------------------------------------------------------
for i = [2:121, 153:158]
    if strcmp(class(data_pd.(variables{i})), 'double')
        idx_nan.PD.(variables{i}) = find(isnan(data_pd.(variables{i})));
    elseif strcmp(class(data_pd.(variables{i})), 'cell')
        idx_nan.PD.(variables{i}) = find(string(data_pd.(variables{i})) == 'NA');
    end
end
clear i

% HC --------------------------------------------------------------------
for i = [2:121, 153:158]
    if strcmp(class(data_hc.(variables{i})), 'double')
        idx_nan.HC.(variables{i}) = find(isnan(data_hc.(variables{i})));
    elseif strcmp(class(data_hc.(variables{i})), 'cell')
        idx_nan.HC.(variables{i}) = find(string(data_hc.(variables{i})) == 'NA');
    end
end
clear i

%%   Handling missing ENROLL_AGE data
% age at the time of the DAT SCAN ---- PD
for i = 1:length(idx_nan.PD.ENROLL_AGE)
    idx_nan_temp = idx_nan.PD.ENROLL_AGE(i);
    datscan_date = char(data_pd.DATSCAN_DATE(idx_nan_temp));
    birth_date = char(data_pd.BIRTHDT(idx_nan_temp));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    data_pd.ENROLL_AGE(idx_nan_temp) = datscan_year - birth_year;
end
idx_nan.PD.ENROLL_AGE = find(isnan(data_pd.ENROLL_AGE));
clear i idx_nan_temp datscan_date datscan_year birth_date birth_year

% age at the time of the DAT SCAN ---- HC
for i = 1:length(idx_nan.HC.ENROLL_AGE)
    idx_nan_temp = idx_nan.HC.ENROLL_AGE(i);
    datscan_date = char(data_hc.DATSCAN_DATE(idx_nan_temp));
    birth_date = char(data_hc.BIRTHDT(idx_nan_temp));
    datscan_year = str2double(string(datscan_date(end-3:end)));
    birth_year = str2double(string(birth_date(end-3:end)));
    data_hc.ENROLL_AGE(idx_nan_temp) = datscan_year - birth_year;
end
idx_nan.HC.ENROLL_AGE = find(isnan(data_hc.ENROLL_AGE));
clear i idx_nan_temp datscan_date datscan_year birth_date birth_year

%% ANALYSIS AVERAGE DEMOGRAPHIC
%SEX
demo.SEX.TOT.num_males = length(data.SEX(strcmpi(data.SEX,'Male')));
demo.SEX.TOT.num_females = length(data.SEX(strcmpi(data.SEX,'Female')));

demo.SEX.PD.num_males = length(data_pd.SEX(strcmpi(data_pd.SEX,'Male')));
demo.SEX.PD.num_females = length(data_pd.SEX(strcmpi(data_pd.SEX,'Female')));

demo.SEX.HC.num_males = length(data_hc.SEX(strcmpi(data_hc.SEX,'Male')));
demo.SEX.HC.num_females = length(data_hc.SEX(strcmpi(data_hc.SEX,'Female')));


% AGE
demo.AGE.TOT.mean = nanmean(data.ENROLL_AGE,"all");
demo.AGE.TOT.min = min(data.ENROLL_AGE);
demo.AGE.TOT.max = max(data.ENROLL_AGE);
demo.AGE.TOT.std = nanstd(data.ENROLL_AGE);


demo.AGE.PD.mean = nanmean(data_pd.ENROLL_AGE,"all");
demo.AGE.PD.min = min(data_pd.ENROLL_AGE);
demo.AGE.PD.max = max(data_pd.ENROLL_AGE);

demo.AGE.HC.mean = nanmean(data_hc.ENROLL_AGE,"all");
demo.AGE.HC.min = min(data_hc.ENROLL_AGE);
demo.AGE.HC.max = max(data_hc.ENROLL_AGE);


%%   Decide what to remove
% Remove from idx_nan empty fields -------- PD
fields_nan_pd = fieldnames(idx_nan.PD);
idx_nan.PD = rmfield(idx_nan.PD, fields_nan_pd(structfun(@isempty, idx_nan.PD)));

% Remove from idx_nan empty fields -------- HC
fields_nan_hc = fieldnames(idx_nan.HC);
idx_nan.HC = rmfield(idx_nan.HC, fields_nan_hc(structfun(@isempty, idx_nan.HC)));

%%    - Remove all PD with NaNs in np test but not np4

idx_fields_tests_no_np4 = [26:54, 61:93];
idx_nan_tests_no_np4 = [];
for i = 1:length(idx_fields_tests_no_np4)
    idx_nan_tests_no_np4 = [idx_nan_tests_no_np4; idx_nan.PD.(fields_nan_pd{idx_fields_tests_no_np4(i)})];
end
idx_nan_tests_no_np4 = unique(idx_nan_tests_no_np4);
clear i idx_fields_tests_no_np4
% Sono 220

% add indeces nan of hand, mcatot, weight, 
idx_nan_all = [idx_nan_tests_no_np4', idx_nan.PD.HANDED', idx_nan.PD.MCATOT', idx_nan.PD.WGTKG'];
idx_nan_all = unique(idx_nan_all);

% Remove patients missing at least one measure of test, excluding np4
data_pd(idx_nan_tests_no_np4,:) = [];

%%    - TRY Remove all HC with NaNs in np test but not np4

idx_fields_tests_no_np4 = [26:54, 60:93];
idx_nan_tests_no_np4 = [];
for i = 1:length(idx_fields_tests_no_np4)
    idx_nan_tests_no_np4 = [idx_nan_tests_no_np4; idx_nan.HC.(fields_nan_pd{idx_fields_tests_no_np4(i)})];
end
idx_nan_tests_no_np4 = unique(idx_nan_tests_no_np4);
clear i idx_fields_tests_no_np4
% Sono 37

% add indeces nan of hand, mcatot, weight, 
idx_nan_all = [idx_nan_tests_no_np4', idx_nan.HC.MCATOT', idx_nan.HC.WGTKG'];
idx_nan_all = unique(idx_nan_all);

% Remove hc missing at least one measure of test, excluding np4
data_hc(idx_nan_tests_no_np4,:) = [];

%%    - Remove patients with low QUALITY RATING
data_pd(find(data_pd.DATSCAN_QUALITY_RATING==3),:) = [];
data_hc(find(data_pd.DATSCAN_QUALITY_RATING==3),:) = [];

%%    - Remove patients with not COMPLETED DATSCAN
data_pd(find(data_pd.DATSCAN==0),:) = [];
data_hc(find(data_pd.DATSCAN==0),:) = [];

%%    - Remove HC patients with abnormal MoCA scores
data_hc(find(data_hc.MCATOT < 26),:) = [];

%% MANAGING OUTLIERS (Box plots)
% remove outliers from both PD nad HC for demographics data

cont_variables = [4, 153, 154]; % age, height, weight
for i = 1:length(cont_variables)
    
    % PD
    [B,TFrm] = rmoutliers(data_pd.(variables{cont_variables(i)}));
    ind_outliers = find(TFrm);
    data_pd(ind_outliers, :) = [];

    % HC
    [B,TFrm] = rmoutliers(data_hc.(variables{cont_variables(i)}));
    ind_outliers = find(TFrm);
    data_hc(ind_outliers, :) = [];
end

clear i B TFrm ind_outliers


%% IMPORTANT VARIABLES EXTRACTION

%%   Dominant hand division
idx_samples.HAND.PD.right = find(strcmp(data_pd.HANDED, 'Right'));
idx_samples.HAND.PD.left = find(strcmp(data_pd.HANDED, 'Left'));
idx_samples.HAND.PD.mixed = find(strcmp(data_pd.HANDED, 'Mixed'));

idx_samples.HAND.HC.right = find(strcmp(data_hc.HANDED, 'Right'));
idx_samples.HAND.HC.left = find(strcmp(data_hc.HANDED, 'Left'));
idx_samples.HAND.HC.mixed = find(strcmp(data_hc.HANDED, 'Mixed'));

%%   Neurophysiological testing and assessment: test in 4 parts
% normalization of np tests scores
np_test_normalized.np1ptot = data_pd.NP1PTOT ./ length(27:33);
np_test_normalized.np1rtot = data_pd.NP1RTOT ./ length(35:40);
np_test_normalized.np2ptot = data_pd.NP2PTOT ./ length(42:54);
np_test_normalized.np3tot = data_pd.NP3TOT ./ length(61:93);
np_test_normalized.np4tot = data_pd.NP4TOT ./ length(98:103);

% % plot per curiosità
% figure, hold on
% plot(np_test_normalized.np1ptot)
% plot(np_test_normalized.np1rtot)
% plot(np_test_normalized.np2ptot)
% plot(np_test_normalized.np3tot)
% plot(np_test_normalized.np4tot)
% axis tight
% xlabel('patients')
% title('normalized NPTEST')
% ylabel('test score')
% legend('1p', '1r', '2p', '3', '4')
% ylim([-0.5 4])

%% ANALYSIS OF DEMOGRAPHICS DATA

%% CONTINUOUS VARIABLES
cont_variables = [4, 153, 154]; % age, height, weight

for i = 1:length(cont_variables)

    % check for gaussianity - PD
    [gaussianity.PD.(variables{cont_variables(i)}).h, gaussianity.PD.(variables{cont_variables(i)}).p] = lillietest(data_pd.(variables{cont_variables(i)}));
    if gaussianity.PD.(variables{cont_variables(i)}).h == 1
        disp(strcat(string(variables{cont_variables(i)}), ' in PD is NOT normally distributed'))
    else
        disp(strcat(string(variables{cont_variables(i)}), ' in PD is normally distributed'))
    end

    % check for gaussianity - HC
    [gaussianity.HC.(variables{cont_variables(i)}).h, gaussianity.HC.(variables{cont_variables(i)}).p] = lillietest(data_hc.(variables{cont_variables(i)}));
    if gaussianity.HC.(variables{cont_variables(i)}).h == 1
        disp(strcat(string(variables{cont_variables(i)}), ' in HC is NOT normally distributed'))
    else
        disp(strcat(string(variables{cont_variables(i)}), ' in HC is normally distributed'))
    end

    % Anova test between PD and HC
    % y = [data_pd.(variables{cont_variables(i)}); data_hc.(variables{cont_variables(i)})]';
    % LAT_groups = [ones(1,size(data_pd.(variables{cont_variables(i)}),1)), 2*ones(1,size(data_hc.(variables{cont_variables(i)}),1))];
    % anova_test.(variables{cont_variables(i)}).p = anova1(y, LAT_groups);
    % 
    % if anova_test.(variables{cont_variables(i)}).p > 0.05
    %     disp(strcat(string(variables{cont_variables(i)}), ': Accepted null hyp (same mean)'))
    % else 
    %     disp(strcat(string(variables{cont_variables(i)}), ': Rejected null hyp (different mean)'))
    % end

end

clear i y LAT_groups

%% DISCRETE VARIABLES
discr_variables = [7, 8, 9, 11]; % ethnicity, sex, hand, familiarity

% Conversion to numerical levels - PD
for i = 1:length(discr_variables)-1 % no need for familiarity
    legend.PD.(variables{discr_variables(i)}) = unique(data_pd.(variables{discr_variables(i)}));
    for j = 1:length(legend.PD.(variables{discr_variables(i)}))
        ind_temp = find(strcmp(data_pd.(variables{discr_variables(i)}), legend.PD.(variables{discr_variables(i)})(j)));
        for k = 1:length(ind_temp)
            data_pd.(variables{discr_variables(i)}){ind_temp(k)} = j;
        end
    end
    data_pd.(variables{discr_variables(i)}) = cellfun(@double, (data_pd.(variables{discr_variables(i)})));
end

% Conversion to numerical levels - HC
for i = 1:length(discr_variables)-1 % no need for familiarity
    legend.HC.(variables{discr_variables(i)}) = unique(data_hc.(variables{discr_variables(i)}));
    for j = 1:length(legend.HC.(variables{discr_variables(i)}))
        ind_temp = find(strcmp(data_hc.(variables{discr_variables(i)}), legend.HC.(variables{discr_variables(i)})(j)));
        for k = 1:length(ind_temp)
            data_hc.(variables{discr_variables(i)}){ind_temp(k)} = j;
        end
    end
    data_hc.(variables{discr_variables(i)}) = cellfun(@double, (data_hc.(variables{discr_variables(i)})));
end

for i = 1:length(discr_variables)
    % Perfom Wilcoxon test
    [wilcoxon.(variables{discr_variables(i)}).p, wilcoxon.(variables{discr_variables(i)}).h] = ranksum(data_pd.(variables{discr_variables(i)}), data_hc.(variables{discr_variables(i)}));
    if wilcoxon.(variables{discr_variables(i)}).h == 1
        disp(strcat(string(variables{discr_variables(i)}), ': Rejected null hyp (different mean)'))
    else
        disp(strcat(string(variables{discr_variables(i)}), ': Accepted null hyp (same mean)'))
    end
end

clear i j k ind_temp

% save legend also for familiarity (ANYFAMPD)
legend.HC.(variables{discr_variables(4)}) = unique(data_hc.(variables{discr_variables(4)}));
legend.HC.(variables{discr_variables(4)}) = legend.HC.(variables{discr_variables(4)})(1:3);

legend.PD.(variables{discr_variables(4)}) = unique(data_pd.(variables{discr_variables(4)}));
legend.PD.(variables{discr_variables(4)}) = legend.PD.(variables{discr_variables(4)})(1:3);


%% LATERALIZATION DATA
% Extraction of lateralization index
% Right = |(right - left)/(right + left)| > 0.20
% Left = |(left - right)/(right + left)| > 0.20
% From: Ipsilateral deficits of dopaminergic neurotransmission in Parkinson s disease

% DATSCAN lateralization PD ------------------------------------------
LATERALIZATION_coeff.CAUDATE.PD = (data_pd.DATSCAN_CAUDATE_R - data_pd.DATSCAN_CAUDATE_L)./(data_pd.DATSCAN_CAUDATE_R + data_pd.DATSCAN_CAUDATE_L);
LATERALIZATION_coeff.PUTAMEN.PD = (data_pd.DATSCAN_PUTAMEN_R - data_pd.DATSCAN_PUTAMEN_L)./(data_pd.DATSCAN_PUTAMEN_R + data_pd.DATSCAN_PUTAMEN_L);
LATERALIZATION_coeff.PUTAMEN_ANT.PD = (data_pd.DATSCAN_PUTAMEN_R_ANT - data_pd.DATSCAN_PUTAMEN_L_ANT)./(data_pd.DATSCAN_PUTAMEN_R_ANT + data_pd.DATSCAN_PUTAMEN_L_ANT);

% DATSCAN lateralization HC ------------------------------------------
LATERALIZATION_coeff.CAUDATE.HC = (data_hc.DATSCAN_CAUDATE_R - data_hc.DATSCAN_CAUDATE_L)./(data_hc.DATSCAN_CAUDATE_R + data_hc.DATSCAN_CAUDATE_L);
LATERALIZATION_coeff.PUTAMEN.HC = (data_hc.DATSCAN_PUTAMEN_R - data_hc.DATSCAN_PUTAMEN_L)./(data_hc.DATSCAN_PUTAMEN_R + data_hc.DATSCAN_PUTAMEN_L);
LATERALIZATION_coeff.PUTAMEN_ANT.HC = (data_hc.DATSCAN_PUTAMEN_R_ANT - data_hc.DATSCAN_PUTAMEN_L_ANT)./(data_hc.DATSCAN_PUTAMEN_R_ANT + data_hc.DATSCAN_PUTAMEN_L_ANT);

% Grafici distribuzioni lateralizzazione?

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.
% With significant difference: 20%

% save indexes of relevant lateralizations - PD
for i = 1:length(ROIs_labels)
    idx_samples.LATERALIZATION.PD.(ROIs_labels(i)).right = find(LATERALIZATION_coeff.(ROIs_labels(i)).PD > 0.2);
    idx_samples.LATERALIZATION.PD.(ROIs_labels(i)).left = find(LATERALIZATION_coeff.(ROIs_labels(i)).PD < -0.2);
    idx_samples.LATERALIZATION.PD.(ROIs_labels(i)).none = find(LATERALIZATION_coeff.(ROIs_labels(i)).PD < 0.2 & LATERALIZATION_coeff.(ROIs_labels(i)).PD > -0.2 );
end

% save indexes of relevant lateralizations - HC
for i = 1:length(ROIs_labels)
    idx_samples.LATERALIZATION.HC.(ROIs_labels(i)).right = find(LATERALIZATION_coeff.(ROIs_labels(i)).HC > 0.2);
    idx_samples.LATERALIZATION.HC.(ROIs_labels(i)).left = find(LATERALIZATION_coeff.(ROIs_labels(i)).HC < -0.2);
    idx_samples.LATERALIZATION.HC.(ROIs_labels(i)).none = find(LATERALIZATION_coeff.(ROIs_labels(i)).HC < 0.2 & LATERALIZATION_coeff.(ROIs_labels(i)).HC > -0.2 );
end

% NOT Absolute value
% DATSCAN lateralization PD ------------------------------------------
NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD = LATERALIZATION_coeff.CAUDATE.PD;
NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD = LATERALIZATION_coeff.PUTAMEN.PD;
NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD = LATERALIZATION_coeff.PUTAMEN_ANT.PD;

% DATSCAN lateralization HC ------------------------------------------
NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC = LATERALIZATION_coeff.CAUDATE.HC;
NOT_ABS.LATERALIZATION_coeff.PUTAMEN.HC = LATERALIZATION_coeff.PUTAMEN.HC;
NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.HC = LATERALIZATION_coeff.PUTAMEN_ANT.HC;


% Absolute value lateralization
% DATSCAN lateralization PD ------------------------------------------
LATERALIZATION_coeff.CAUDATE.PD = abs(LATERALIZATION_coeff.CAUDATE.PD);
LATERALIZATION_coeff.PUTAMEN.PD = abs(LATERALIZATION_coeff.PUTAMEN.PD);
LATERALIZATION_coeff.PUTAMEN_ANT.PD = abs(LATERALIZATION_coeff.PUTAMEN_ANT.PD);

% DATSCAN lateralization HC ------------------------------------------
LATERALIZATION_coeff.CAUDATE.HC = abs(LATERALIZATION_coeff.CAUDATE.HC);
LATERALIZATION_coeff.PUTAMEN.HC = abs(LATERALIZATION_coeff.PUTAMEN.HC);
LATERALIZATION_coeff.PUTAMEN_ANT.HC = abs(LATERALIZATION_coeff.PUTAMEN_ANT.HC);

%% VISUALISING LATERALIZATION
% difference_caudate_hc = mean(LATERALIZATION_coeff.CAUDATE.HC);
% difference_caudate_pd = mean(LATERALIZATION_coeff.CAUDATE.PD);
% 
% difference_putamen_hc = mean(LATERALIZATION_coeff.PUTAMEN.HC);
% difference_putamen_pd = mean(LATERALIZATION_coeff.PUTAMEN.PD);
% 
% difference_putamen_ant_hc = mean(LATERALIZATION_coeff.PUTAMEN_ANT.HC);
% difference_putamen_ant_pd = mean(LATERALIZATION_coeff.PUTAMEN_ANT.PD);
% 
% difference_hc = [NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.HC,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.HC];
% difference_pd = [NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD];
% 
% figure
% histogram(difference_hc,'Normalization','countdensity') 
% hold on
% % ylim([-0.8, 0.8])
% % pd
% % boxchart(difference_pd, 'MarkerStyle', 'none', 'BoxFaceColor', 'r')
% histogram(difference_pd,'Normalization','countdensity')
% % visualization
% % view(-90,90)
% % xlabel('ROIs')
% % ylabel('DAT difference')
% % ylim([-0.8, 0.8])
% title('Difference DAT_{left} - DAT_{right} for HC and PD')
% % xticklabels({'Caudate','Putamen','Putamen Anterior'})
% % yline(0)
% legend('HC','PD')
% 
% 
difference_caudate_hc = mean(LATERALIZATION_coeff.CAUDATE.HC);
difference_caudate_pd = mean(LATERALIZATION_coeff.CAUDATE.PD);

difference_putamen_hc = mean(LATERALIZATION_coeff.PUTAMEN.HC);
difference_putamen_pd = mean(LATERALIZATION_coeff.PUTAMEN.PD);

difference_putamen_ant_hc = mean(LATERALIZATION_coeff.PUTAMEN_ANT.HC);
difference_putamen_ant_pd = mean(LATERALIZATION_coeff.PUTAMEN_ANT.PD);

difference_hc = [NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.HC,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.HC];
difference_pd = [NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD];

figure
histogram(difference_hc,'Normalization','countdensity') 
hold on
% ylim([-0.8, 0.8])
% pd
% boxchart(difference_pd, 'MarkerStyle', 'none', 'BoxFaceColor', 'r')
histogram(difference_pd,'Normalization','countdensity')
% visualization
% view(-90,90)
% xlabel('ROIs')
% ylabel('DAT difference')
% ylim([-0.8, 0.8])
title('Difference DAT_{left} - DAT_{right} for HC and PD')
% xticklabels({'Caudate','Putamen','Putamen Anterior'})
% yline(0)
% legend('HC','PD')



%% STATISTICAL ANALYSIS 
%% One-way Anova
% Control gaussianity
rois_names = fieldnames(LATERALIZATION_coeff);
for i = 1:length(rois_names)
    rois_names_diag = fieldnames(LATERALIZATION_coeff.(rois_names{i,1}));
    for j = 1:length(rois_names_diag)
        [gaussianity.(rois_names_diag{j}).(rois_names{i}).h, gaussianity.(rois_names_diag{j}).(rois_names{i}).p] = lillietest(LATERALIZATION_coeff.(rois_names{i}).(rois_names_diag{j}));
        if gaussianity.(rois_names_diag{j}).(rois_names{i}).h == 1
            disp([rois_names{i}; rois_names_diag{j}; "not gaussian"])
        end
    end
end
clear rois_names rois_names_diag i j

% analysis with absolute value
% Caudate
y = [LATERALIZATION_coeff.CAUDATE.PD' LATERALIZATION_coeff.CAUDATE.HC']';
group1 = [repmat("PD_caudate",length(LATERALIZATION_coeff.CAUDATE.PD),1);repmat("HC_caudate",length(LATERALIZATION_coeff.CAUDATE.HC),1)];
[p,tbl,stats]  = anova1(y,group1);
if p > 0.05
    disp("There isn't significant difference between the means of the DATSCAN of HC and PD in the CAUDATE region")
else
    disp("There is significant difference between the means of the DATSCAN of HC and PD in the CAUDATE region")
end

% Putamen
y = [LATERALIZATION_coeff.PUTAMEN.PD' LATERALIZATION_coeff.PUTAMEN.HC']';
group2 = [repmat("PD_putamen",length(LATERALIZATION_coeff.PUTAMEN.PD),1);repmat("HC_putamen",length(LATERALIZATION_coeff.PUTAMEN.HC),1)];
[p,tbl,stats]  = anova1(y,group2);
if p > 0.05
    disp("There isn't significant difference between the means of the DATSCAN of HC and PD in the PUTAMEN region")
else
    disp("There is significant difference between the means of the DATSCAN of HC and PD in the PUTAMEN region")
end

% Putamen ANT
y = [LATERALIZATION_coeff.PUTAMEN_ANT.PD' LATERALIZATION_coeff.PUTAMEN_ANT.HC']';
group3 = [repmat("PD_putamen_ant",length(LATERALIZATION_coeff.PUTAMEN_ANT.PD),1);repmat("HC_putamen_ant",length(LATERALIZATION_coeff.PUTAMEN_ANT.HC),1)];
[p,tbl,stats]  = anova1(y,group3);
if p > 0.05
    disp("There isn't significant difference between the means of the DATSCAN of HC and PD in the PUTAMEN ANT region")
else
    disp("There is significant difference between the means of the DATSCAN of HC and PD in the PUTAMEN ANT region")
end

% with all the ROIs
y = [LATERALIZATION_coeff.CAUDATE.PD' LATERALIZATION_coeff.CAUDATE.HC' LATERALIZATION_coeff.PUTAMEN.PD' LATERALIZATION_coeff.PUTAMEN.HC' LATERALIZATION_coeff.PUTAMEN_ANT.PD' LATERALIZATION_coeff.PUTAMEN_ANT.HC' ]';
group4 = [group1; group2; group3];
[p,tbl,stats]  = anova1(y,group4);
title('Lateralization coefficients for each ROI - PD & HC')
ylabel('Lateralization coefficient')

clear y p tbl stats group4 group3 group2 group1 



%% COVARIATES ANALYSIS
% Comparison of lateralization coefficient only within SEX 
% The other variables don't have enough samples

cohorts = ["HC", "PD"];
for cohort = cohorts
    for region = ROIs_labels
        for i = 1:length(discr_variables)-1

            % groups LAT
            for j = 1:length(legend.(cohort).(variables{discr_variables(i)}))
                LAT_groups.(cohort).(region).(legend.(cohort).(variables{discr_variables(i)}){j}) = LATERALIZATION_coeff.(region).(cohort)(find(data_hc.(variables{discr_variables(i)}) == j));

                % not consider groups with less than 4 elements
                if length(LAT_groups.(cohort).(region).(legend.(cohort).(variables{discr_variables(i)}){j})) > 4
                    % check for gaussianity
                    [gaussianity.(cohort).LATERALIZATION.(region).(variables{discr_variables(i)}).(legend.(cohort).(variables{discr_variables(i)}){j}).h, ...
                        gaussianity.(cohort).LATERALIZATION.(region).(variables{discr_variables(i)}).(legend.(cohort).(variables{discr_variables(i)}){j}).p] = ...
                        lillietest(LAT_groups.(cohort).(region).(legend.(cohort).(variables{discr_variables(i)}){j}));
                end
            end
        end

        % Comparison region lateralization female - male in HC
        % Anova test
        y = [LAT_groups.(cohort).(region).Female; LAT_groups.HC.(region).Male]';
        groups_anova = [ones(1,length(LAT_groups.(cohort).(region).Female)), 2*ones(1,length(LAT_groups.(cohort).(region).Male))];
        anova_test.LATERALIZATION.(region).(cohort).(variables{discr_variables(2)}).p = anova1(y, groups_anova);

        if anova_test.LATERALIZATION.(region).(cohort).(variables{discr_variables(2)}).p > 0.05
            disp(strcat(cohort, string(variables{discr_variables(2)}), ': LAT', region, ' Accepted null hyp (same mean)'))
        else 
            disp(strcat(cohort, string(variables{discr_variables(2)}), ': LAT', region, ' Rejected null hyp (different mean)'))
        end

        wilcoxon_test.LATERALIZATION.(region).(cohort).(variables{discr_variables(2)}).p = ranksum(LAT_groups.(cohort).(region).Female, LAT_groups.HC.(region).Male);

        if wilcoxon_test.LATERALIZATION.(region).(cohort).(variables{discr_variables(2)}).p > 0.05
            disp(strcat(cohort, string(variables{discr_variables(2)}), ': LAT ', region, ' Accepted null hyp (same mean)'))
        else 
            disp(strcat(cohort,string(variables{discr_variables(2)}), ': LAT ', region, ' Rejected null hyp (different mean)'))
        end

        clear y groups_anova

    end
    clear i j region
end


% % Difference in PUTAMEN and PUTAMEN ANT lateralization between male and
% % female
% figure
% boxplot(LAT_groups.HC.PUTAMEN.Male, 'Colors','r')
% hold on
% boxplot(LAT_groups.HC.PUTAMEN.Female,'Colors','b')
% 
% figure
% boxplot(LAT_groups.HC.PUTAMEN_ANT.Male, 'Colors','r')
% hold on
% boxplot(LAT_groups.HC.PUTAMEN_ANT.Female,'Colors','b')


%% LINEAR REGRESSION of variables of interest



%% CORRELATION MATRIX of variables of interest
% fare scatter
%%   - HC - indirect data
idx_covariates = [4, 153, 154]; 
covariates_hc = data_hc(:, idx_covariates); % age, height, weight
covariates_hc.('CAUDATE') = LATERALIZATION_coeff.CAUDATE.HC;
covariates_hc.('PUTAMEN') = LATERALIZATION_coeff.PUTAMEN.HC;
covariates_hc.('PUTAMEN_ANT') = LATERALIZATION_coeff.PUTAMEN_ANT.HC;

[correlation.indirect.HC.p, correlation.indirect.HC.R] = corrcoef(table2array(covariates_hc), 'Rows', 'pairwise');

figure
imagesc(correlation.indirect.HC.R)
colormap parula
colorbar
xticks(1:width(covariates_hc))
xticklabels(["ENROLL\_AGE", "HTCM", "WGTKG", "CAUDATE", "PUTAMEN", "PUTAMEN\_ANT"])
yticks(1:width(covariates_hc))
yticklabels(["ENROLL\_AGE", "HTCM", "WGTKG", "CAUDATE", "PUTAMEN", "PUTAMEN\_ANT"])
title("Correlation between lateralization and age, weight, height")
ax = gca;
ax.FontSize = 9;
axis equal
for i=1:3
    covariates_to_save_hc.(ROIs_labels(i)) = [];
end

for i =1:size(correlation.indirect.HC.R,1)
    for j = 1:size(correlation.indirect.HC.R,2)
        if correlation.indirect.HC.p(i,j) < 0.05 && correlation.indirect.HC.R(i,j) > 0.5
                if (j <= (size(correlation.indirect.HC.R,1) -3)) && (i > (size(correlation.indirect.HC.R,1) -3)) 
                    disp([covariates_hc.Properties.VariableNames{i}, ' correlated with ', covariates_hc.Properties.VariableNames{j}])
                    covariates_to_save_hc.(ROIs_labels(i-3)) = [covariates_to_save_hc.(ROIs_labels(i-3)),  convertCharsToStrings(covariates_hc.Properties.VariableNames{j})];
                end
       end
    end
end

for i=1:3
    covariates_to_save_hc.(ROIs_labels(i))  = unique(covariates_to_save_hc.(ROIs_labels(i)));
end
%%   - HC - symptoms data
idx_symptoms = [34, 55, 94, 27:54, 61:93,158]; 
symptoms_hc = data_hc(:, idx_symptoms); % age, height, weight, mcatot
symptoms_hc.('CAUDATE') = LATERALIZATION_coeff.CAUDATE.HC;
symptoms_hc.('PUTAMEN') = LATERALIZATION_coeff.PUTAMEN.HC;
symptoms_hc.('PUTAMEN_ANT') = LATERALIZATION_coeff.PUTAMEN_ANT.HC;

[correlation.symptoms.HC.p, correlation.symptoms.HC.R] = corrcoef(table2array(symptoms_hc), 'Rows', 'pairwise');

figure
imagesc(correlation.symptoms.HC.R)
colormap parula
colorbar
xticks(1:width(symptoms_hc))
xticklabels(symptoms_hc.Properties.VariableNames)
yticks(1:width(symptoms_hc))
yticklabels(symptoms_hc.Properties.VariableNames)
title("Correlation between lateralization and symptoms - HC")
ax = gca;
ax.FontSize = 6;
axis equal

for i=1:3
    symptoms_to_save_hc.(ROIs_labels(i)) = [];
end


for i =1:size(correlation.symptoms.HC.R,1)
    for j = 1:size(correlation.symptoms.HC.R,2)
        if correlation.symptoms.HC.p(i,j) < 0.05 && correlation.symptoms.HC.R(i,j) > 0.5
                if (j <= (size(correlation.symptoms.HC.R,1) -3)) && (i > (size(correlation.symptoms.HC.R,1) -3)) 
                    disp([symptoms_hc.Properties.VariableNames{i}, ' correlated with ', symptoms_hc.Properties.VariableNames{j}])
                    symptoms_to_save_hc.(ROIs_labels(size(correlation.symptoms.HC.R,1) -i+1)) = [symptoms_to_save_hc.(ROIs_labels(size(correlation.symptoms.HC.R,1) -i+1)),  convertCharsToStrings(symptoms_hc.Properties.VariableNames{j})];
                end
       end
    end
end

for i=1:3
    symptoms_to_save_hc.(ROIs_labels(i))  = unique(symptoms_to_save_hc.(ROIs_labels(i)));
end



%%  - PD - indirect data
idx_tot = [idx_covariates,idx_symptoms];
covariates_pd = data_pd(:, idx_covariates); % age, weight, height,  np test + mcatot
covariates_pd.('CAUDATE') = LATERALIZATION_coeff.CAUDATE.PD;
covariates_pd.('PUTAMEN') = LATERALIZATION_coeff.PUTAMEN.PD;
covariates_pd.('PUTAMEN_ANT') = LATERALIZATION_coeff.PUTAMEN_ANT.PD;

[correlation.indirect.PD.p, correlation.indirect.PD.R] = corrcoef(table2array(covariates_pd), 'Rows', 'pairwise');

figure
imagesc(correlation.indirect.PD.R)
colormap parula
colorbar
xticks(1:width(covariates_pd))
xticklabels(covariates_pd.Properties.VariableNames)
yticks(1:width(covariates_pd))
yticklabels(covariates_pd.Properties.VariableNames)
title("Correlation between lateralization and age, weight, height, MCATOT - PD")
ax = gca;
ax.FontSize = 6;
axis equal

for i=1:3
    covariates_to_save_pd.(ROIs_labels(i)) = [];
end

for i =1:size(correlation.indirect.PD.R,1)
    for j = 1:size(correlation.indirect.PD.R,2)
        if correlation.indirect.PD.p(i,j) < 0.05 && correlation.indirect.PD.R(i,j) > 0.5
                if (j < (size(correlation.indirect.PD.R,1) -3)) && (i > (size(correlation.indirect.PD.R,1) -3)) 
                    disp([covariates_pd.Properties.VariableNames{i}, ' correlated with ', covariates_pd.Properties.VariableNames{j}])
                    covariates_to_save_pd.(ROIs_labels(i-3)) = [covariates_to_save_pd.(ROIs_labels(i-3)),  convertCharsToStrings(covariates_pd.Properties.VariableNames{j})];
                end
       end
    end
end


for i=1:3
    covariates_to_save_pd.(ROIs_labels(i))  = unique(covariates_to_save_pd.(ROIs_labels(i)));
end


%%   - PD - symptoms data
idx_symptoms = [34, 55, 94, 27:54, 61:93]; 
symptoms_pd = data_pd(:, idx_symptoms); % age, height, weight, mcatot
symptoms_pd.('CAUDATE') = LATERALIZATION_coeff.CAUDATE.PD;
symptoms_pd.('PUTAMEN') = LATERALIZATION_coeff.PUTAMEN.PD;
symptoms_pd.('PUTAMEN_ANT') = LATERALIZATION_coeff.PUTAMEN_ANT.PD;

[correlation.symptoms.PD.p, correlation.symptoms.PD.R] = corrcoef(table2array(symptoms_pd), 'Rows', 'pairwise');

figure
imagesc(correlation.symptoms.PD.R)
colormap parula
colorbar
xticks(1:width(symptoms_pd))
xticklabels(symptoms_pd.Properties.VariableNames)
yticks(1:width(symptoms_pd))
yticklabels(symptoms_pd.Properties.VariableNames)
t = title("Correlation between lateralization and symptoms - PD");
ax = gca;
ax.FontSize = 5;
axis equal
t.FontSize = 15;

for i=1:3
    symptoms_to_save_pd.(ROIs_labels(i)) = [];
end

for i =1:size(correlation.symptoms.PD.R,1)
    for j = 1:size(correlation.symptoms.PD.R,2)
        if correlation.symptoms.PD.p(i,j) < 0.05 && correlation.symptoms.PD.R(i,j) > 0.5
                if (j <= (size(correlation.symptoms.PD.R,1) -3)) && (i > (size(correlation.symptoms.PD.R,1) -3)) 
                    disp([symptoms_pd.Properties.VariableNames{i}, ' correlated with ', symptoms_pd.Properties.VariableNames{j}])
                    symptoms_to_save_pd.(ROIs_labels(size(correlation.symptoms.HC.R,1) -i)) = [symptoms_to_save_pd.(ROIs_labels(size(correlation.symptoms.HC.R,1) -i)),  convertCharsToStrings(symptoms_pd.Properties.VariableNames{j})];
                end
       end
    end
end


for i=1:3
    symptoms_to_save_pd.(ROIs_labels(i))  = unique(symptoms_to_save_pd.(ROIs_labels(i)));
end

%% LINEAR REGRESSION
%% - HC

for i=1:3
    symptoms_hc_array.(ROIs_labels(i)) = table2array(symptoms_hc(:,contains(symptoms_hc.Properties.VariableNames,symptoms_to_save_hc.(ROIs_labels(i)))));
    try
        covariates_hc_array.(ROIs_labels(i)) = table2array(covariates_hc(:,contains(covariates_hc.Properties.VariableNames,covariates_to_save_hc.(ROIs_labels(i)))));
    catch
        covariates_hc_array.(ROIs_labels(i)) = [];
    end
    regressors_hc_array.(ROIs_labels(i)) = [covariates_hc_array.(ROIs_labels(i)),symptoms_hc_array.(ROIs_labels(i))];

end

figure(25)
set(gcf, 'Position', get(0, 'Screensize'));

subplot(131)
model_caud_hc = fitlm(covariates_hc_array.CAUDATE,NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC);
plot(model_caud_hc)
% xlim([-0.05 0.05])
% ylim([-0.3 0.2])
xlabel('Covariates')
ylabel('Lateralization coefficient')
title('Caudate linear fit - HC')

subplot(132)
model_put_hc = fitlm(covariates_hc_array.PUTAMEN,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.HC);
plot(model_put_hc)
% xlim([-0.05 0.05])
% ylim([-0.3 0.2])
xlabel('Covariates')
ylabel('Lateralization coefficient')
title('Putamen linear fit - HC')

subplot(133)
model_put_ant_hc = fitlm(covariates_hc_array.PUTAMEN_ANT,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.HC);
plot(model_put_ant_hc)
% xlim([-0.05 0.05])
% ylim([-0.3 0.2])
xlabel('Covariates')
ylabel('Lateralization coefficient')
title('Putamen Anterior linear fit - HC')
sgtitle('Lateralization coefficient - Covariates FIT HC')


saveas(figure(25), "fit_covariates_hc.png", "png")

%% ------ Statistics
anova_caud_hc = anova(model_caud_hc,'component');
anova_put_hc = anova(model_put_hc,'summary');
anova_put_ant_hc =  anova(model_put_ant_hc,'summary');

%% - PD

for i=1:3
    symptoms_pd_array.(ROIs_labels(i)) = table2array(symptoms_pd(:,contains(symptoms_pd.Properties.VariableNames,symptoms_to_save_pd.(ROIs_labels(i)))));
    try
        covariates_pd_array.(ROIs_labels(i)) = table2array(covariates_pd(:,contains(covariates_pd.Properties.VariableNames,covariates_to_save_hc.(ROIs_labels(i)))));
    catch
        covariates_pd_array.(ROIs_labels(i)) = [];
    end

    regressors_pd_array.(ROIs_labels(i)) = [covariates_pd_array.(ROIs_labels(i)),symptoms_pd_array.(ROIs_labels(i))];

end

a = regressors_pd_array.CAUDATE;
a(:,2) = [];
b = regressors_pd_array.CAUDATE;
b(:,3) = [];
c = [regressors_pd_array.CAUDATE(:,2),regressors_pd_array.CAUDATE(:,4:end)];
d = [regressors_pd_array.CAUDATE(:,1),regressors_pd_array.CAUDATE(:,4:end)];



permutazioni_caud = {regressors_pd_array.CAUDATE; regressors_pd_array.CAUDATE(:,2:end); regressors_pd_array.CAUDATE(:,3:end);regressors_pd_array.CAUDATE(:,4:end);a ; b;c;d};
names_permutazioni_caud = {'ALL';'HEIGHT & WEIGHT';'WEIGHT';'NOTHING';'AGE & WEIGHT';'AGE & HEIGHT';'AGE';'HEIGHT'};
R_squared_caud = [];
model_caud_pd = {};

clear a b c d

permutazioni_put = {regressors_pd_array.PUTAMEN;regressors_pd_array.PUTAMEN(:,2:end)};
names_permutazioni_put = {'WEIGHT';'NOTHING'};
R_squared_put = [];
model_put_pd= {};

a = regressors_pd_array.PUTAMEN_ANT;
a(:,2) = [];
permutazioni_put_ant = {regressors_pd_array.PUTAMEN_ANT; regressors_pd_array.PUTAMEN_ANT(:,2:end); regressors_pd_array.PUTAMEN_ANT(:,3:end); a};
names_permutazioni_put_ant = {'ALL';'WEIGHT';'NOTHING';'HEIGHT'};
R_squared_put_ant = [];
model_put_ant_pd = {};

clear a



figure(26)
set(gcf, 'Position', get(0, 'Screensize'));

subplot(131)
for i=1:8
    model_caud_pd{i} = fitlm(permutazioni_caud{i},NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
    R_squared_caud(i) = model_caud_pd{i}.Rsquared.Adjusted;
end
[max_R_squared_caud,idx_max_R_squared_caud] = max(R_squared_caud);
plot(model_caud_pd{idx_max_R_squared_caud})
covariate_fit.CAUDATE.ALL = names_permutazioni_caud{idx_max_R_squared_caud};
xlim([-4 6])
ylim([-0.6 0.8])
xlabel('Covariates and symptoms')
ylabel('Lateralization coefficient')
title('Caudate linear fit - PD')

subplot(132)
for i=1:2
    model_put_pd{i} = fitlm(permutazioni_put{i},NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD);
    R_squared_put(i) = model_put_pd{i}.Rsquared.Adjusted;
end
[max_R_squared_put,idx_max_R_squared_put] = max(R_squared_put);

plot(model_put_pd{idx_max_R_squared_put})
covariate_fit.PUTAMEN.ALL = names_permutazioni_put{idx_max_R_squared_put};


xlim([-4 6])
ylim([-0.6 0.8])
xlabel('Covariates and symptoms')
ylabel('Lateralization coefficient')
title('Putamen linear fit - PD')

subplot(133)
for i=1:4
    model_put_ant_pd{i} = fitlm(permutazioni_put_ant{i},NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD);
    R_squared_put_ant(i) = model_put_ant_pd{i}.Rsquared.Adjusted;
end
[max_R_squared_put_ant,idx_max_R_squared_put_ant] = max(R_squared_put_ant);

plot(model_put_ant_pd{idx_max_R_squared_put_ant})
covariate_fit.PUTAMEN_ANT.ALL = names_permutazioni_put_ant{idx_max_R_squared_put_ant};


xlim([-4 6])
ylim([-0.6 0.8])
xlabel('Covariates and symptoms')
ylabel('Lateralization coefficient')
title('Putamen Anterior linear fit - PD')
sgtitle('Lateralization coefficient - Regressors FIT PD')


saveas(figure(26), "fit_covariates_pd.png", "png")

%% ----- Statistics
anova_caud_pd = anova(model_caud_pd{idx_max_R_squared_caud},'summary');
anova_put_pd = anova(model_put_pd{idx_max_R_squared_put},'summary');
anova_put_ant_pd = anova(model_put_ant_pd{idx_max_R_squared_put},'summary');


%% SEPARATIONS Male - Female
idx_male_population = find(data_pd.SEX== 2);
idx_female_population = find(data_pd.SEX== 1);

% male_lateralization_coeff.CAUDATE = NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD(idx_male_population,:);
% female_lateralization_coeff.CAUDATE = NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD(idx_female_population,:);
% 
% male_lateralization_coeff.PUTAMEN = NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD(idx_male_population,:);
% female_lateralization_coeff.PUTAMEN = NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD(idx_female_population,:);

male_lateralization_coeff.PUTAMEN_ANT = NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD(idx_male_population,:);
female_lateralization_coeff.PUTAMEN_ANT = NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD(idx_female_population,:);

for i=1:3
    regressors_male_array.(ROIs_labels(i)) = regressors_pd_array.(ROIs_labels(i))(idx_male_population,:);
    regressors_female_array.(ROIs_labels(i)) = regressors_pd_array.(ROIs_labels(i))(idx_female_population,:);
end

% Permutazioni Male
% a = regressors_male_array.CAUDATE;
% a(:,2) = [];
% b = regressors_male_array.CAUDATE;
% b(:,3) = [];
% c = [regressors_male_array.CAUDATE(:,2),regressors_male_array.CAUDATE(:,4:end)];
% d = [regressors_male_array.CAUDATE(:,1),regressors_male_array.CAUDATE(:,4:end)];
% 
% permutazioni_caud_male = {regressors_male_array.CAUDATE; regressors_male_array.CAUDATE(:,2:end); regressors_male_array.CAUDATE(:,3:end);regressors_male_array.CAUDATE(:,4:end);a ; b;c;d};
% R_squared_caud_male = [];
% model_caud_male = {};
% 
% clear a b c d
% 
% permutazioni_put_male = {regressors_male_array.PUTAMEN;regressors_male_array.PUTAMEN(:,2:end)};
% R_squared_put_male = [];
% model_put_male= {};
% 
a = regressors_male_array.PUTAMEN_ANT;
a(:,2) = [];
permutazioni_put_ant_male = {regressors_male_array.PUTAMEN_ANT; regressors_male_array.PUTAMEN_ANT(:,2:end); regressors_male_array.PUTAMEN_ANT(:,3:end); a};
R_squared_put_ant_male = [];
model_put_ant_male = {};

clear a


% Male fit
% for i=1:8
%     model_caud_male{i} = fitlm(permutazioni_caud_male{i},male_lateralization_coeff.CAUDATE);
%     R_squared_caud_male(i) = model_caud_male{i}.Rsquared.Adjusted;
% end
% [max_R_squared_caud_male,idx_max_R_squared_caud_male] = max(R_squared_caud_male);
% covariate_fit.CAUDATE.MALE = names_permutazioni_caud{idx_max_R_squared_caud_male};
% 
% for i=1:2
%     model_put_male{i} = fitlm(permutazioni_put_male{i},male_lateralization_coeff.PUTAMEN);
%     R_squared_put_male(i) = model_put_male{i}.Rsquared.Adjusted;
% end
% [max_R_squared_put_male,idx_max_R_squared_put_male] = max(R_squared_put_male);
% covariate_fit.PUTAMEN.MALE = names_permutazioni_put{idx_max_R_squared_put_male};
% 

for i=1:4
    model_put_ant_male{i} = fitlm(permutazioni_put_ant_male{i},male_lateralization_coeff.PUTAMEN_ANT);
    R_squared_put_ant_male(i) = model_put_ant_male{i}.Rsquared.Adjusted;
end
[max_R_squared_put_ant_male,idx_max_R_squared_put_ant_male] = max(R_squared_put_ant_male);
covariate_fit.PUTAMEN_ANT.MALE = names_permutazioni_put_ant{idx_max_R_squared_put_ant_male};


figure(27)
set(gcf, 'Position', get(0, 'Screensize'));

% subplot(131)
% plot(model_caud_male{idx_max_R_squared_caud_male})
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% t = title('Fit Lateralization CAUDATE - Regressors [Male population]');
% set(t, 'FontSize', 9)
% subplot(132)
% plot(model_put_male{idx_max_R_squared_put_male})
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% t = title('Fit Lateralization PUTAMEN - Regressors [Male population]');
% set(t, 'FontSize', 9)
% subplot(133)




% Permutazioni Female
% a = regressors_female_array.CAUDATE;
% a(:,2) = [];
% b = regressors_female_array.CAUDATE;
% b(:,3) = [];
% c = [regressors_female_array.CAUDATE(:,2),regressors_female_array.CAUDATE(:,4:end)];
% d = [regressors_female_array.CAUDATE(:,1),regressors_female_array.CAUDATE(:,4:end)];
% 
% permutazioni_caud_female = {regressors_female_array.CAUDATE; regressors_female_array.CAUDATE(:,2:end); regressors_female_array.CAUDATE(:,3:end);regressors_female_array.CAUDATE(:,4:end);a ; b;c;d};
% R_squared_caud_female = [];
% model_caud_female = {};
% 
% clear a b c d
% 
% permutazioni_put_female = {regressors_female_array.PUTAMEN;regressors_female_array.PUTAMEN(:,2:end)};
% R_squared_put_female = [];
% model_put_female= {};
% 
a = regressors_female_array.PUTAMEN_ANT;
a(:,2) = [];
permutazioni_put_ant_female = {regressors_female_array.PUTAMEN_ANT; regressors_female_array.PUTAMEN_ANT(:,2:end); regressors_female_array.PUTAMEN_ANT(:,3:end); a};
R_squared_put_ant_female = [];
model_put_ant_female = {};

clear a

% Female fit

% for i=1:8
%     model_caud_female{i} = fitlm(permutazioni_caud_female{i},female_lateralization_coeff.CAUDATE);
%     R_squared_caud_female(i) = model_caud_female{i}.Rsquared.Adjusted;
% end
% [max_R_squared_caud_female,idx_max_R_squared_caud_female] = max(R_squared_caud_female);
% covariate_fit.CAUDATE.FEMALE = names_permutazioni_caud{idx_max_R_squared_caud_female};
% 
% for i=1:2
%     model_put_female{i} = fitlm(permutazioni_put_female{i},female_lateralization_coeff.PUTAMEN);
%     R_squared_put_female(i) = model_put_female{i}.Rsquared.Adjusted;
% end
% [max_R_squared_put_female,idx_max_R_squared_put_female] = max(R_squared_put_female);
% covariate_fit.PUTAMEN.FEMALE = names_permutazioni_put{idx_max_R_squared_put_female};

for i=1:4
    model_put_ant_female{i} = fitlm(permutazioni_put_ant_female{i},female_lateralization_coeff.PUTAMEN_ANT);
    R_squared_put_ant_female(i) = model_put_ant_female{i}.Rsquared.Adjusted;
end
[max_R_squared_put_ant_female,idx_max_R_squared_put_ant_female] = max(R_squared_put_ant_female);
covariate_fit.PUTAMEN_ANT.FEMALE = names_permutazioni_put_ant{idx_max_R_squared_put_ant_female};

% figure(28)
% set(gcf, 'Position', get(0, 'Screensize'));

% subplot(131)
% plot(model_caud_female{idx_max_R_squared_caud_female})
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% t = title('Fit Lateralization CAUDATE - Regressors [Female population]');
% set(t, 'FontSize', 9)
% subplot(132)
% plot(model_put_female{idx_max_R_squared_put_female})
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% t = title('Fit Lateralization PUTAMEN - Regressors [Female population]');
% set(t, 'FontSize', 9)
% subplot(133)
figure(27)
subplot(121)
plot(model_put_ant_male{idx_max_R_squared_put_ant_male})
xlabel('Covariates and symptoms')
ylabel('Lateralization coefficient')
ylim([-0.4 0.5])
xlim([-2.5 7])
t = title('Fit Lateralization PUTAMEN ANTERIOR - Regressors [Male population]');
set(t, 'FontSize', 9)

% sgtitle('Male population FIT')
subplot(122)
plot(model_put_ant_female{idx_max_R_squared_put_ant_female})
xlabel('Covariates and symptoms')
ylabel('Lateralization coefficient')
ylim([-0.4 0.5])
xlim([-2.5 7])
t = title('Fit Lateralization PUTAMEN ANTERIOR - Regressors [Female population]');
set(t, 'FontSize', 9)

sgtitle('Male and female population FIT')
saveas(figure(27), "fit_covariates_male_female.png", "png")
% clear t



%% SIMPTOMS LATERALITY INDEX

% rigidity_upper = (data_pd.NP3RIGRU -data_pd.NP3RIGLU)./(data_pd.NP3RIGRU + data_pd.NP3RIGLU);  
% rigidity_lower = (data_pd.NP3RIGRL -data_pd.NP3RIGLL)./(data_pd.NP3RIGRL + data_pd.NP3RIGLL);
% tap_hand = (data_pd.NP3FTAPR - data_pd.NP3FTAPL)./(data_pd.NP3FTAPR + data_pd.NP3FTAPL);
% move_hand = (data_pd.NP3HMOVL - data_pd.NP3HMOVR)./(data_pd.NP3HMOVL + data_pd.NP3HMOVR); 
% pu_hand = (data_pd.NP3PRSPR - data_pd.NP3PRSPL)./(data_pd.NP3PRSPR + data_pd.NP3PRSPL);  
% tap_foot = (data_pd.NP3TTAPR - data_pd.NP3TTAPL)./(data_pd.NP3TTAPR + data_pd.NP3TTAPL);
% tap_leg = (data_pd.NP3LGAGR - data_pd.NP3LGAGL)./(data_pd.NP3LGAGR + data_pd.NP3LGAGL);
% post_trem_leg =  (data_pd.NP3PTRMR - data_pd.NP3PTRML)./((data_pd.NP3PTRMR + data_pd.NP3PTRML));
% kin_trem_hand = (data_pd.NP3KTRMR-data_pd.NP3KTRML)./(data_pd.NP3KTRMR + data_pd.NP3KTRML);
% rest_trem_up = (data_pd.NP3RTARU - data_pd.NP3RTALU)./(data_pd.NP3RTARU + data_pd.NP3RTALU);
% rest_trem_low = (data_pd.NP3RTARL - data_pd.NP3RTALL)./(data_pd.NP3RTARL + data_pd.NP3RTALL);
% 
% 
% sintomi_lat = [rigidity_upper';rigidity_lower';tap_hand';move_hand';pu_hand';tap_foot';tap_leg';post_trem_leg';kin_trem_hand';rest_trem_up';rest_trem_low']';
% 
% model_sintomi_lat_caud = fitlm(sintomi_lat,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
% model_sintomi_lat_put = fitlm(sintomi_lat,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD);
% model_sintomi_lat_put_ant = fitlm(sintomi_lat,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD);
% 
% figure(29)
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% subplot(131)
% plot(model_sintomi_lat_caud)
% xlabel('Asymmetry index symptoms')
% ylabel('Lateralization index')
% title('Caudate linear fit - PD')
% 
% subplot(132)
% plot(model_sintomi_lat_put)
% xlabel('Asymmetry index symptoms')
% ylabel('Lateralization index')
% title('Putamen linear fit - PD')
% 
% subplot(133)
% plot(model_sintomi_lat_put_ant)
% xlabel('Asymmetry index symptoms')
% ylabel('Lateralization index')
% title('Putamen Anterior linear fit - PD')
% 
% saveas(figure(27), "fit asimmetry index.png", "png")
% 





% %% AGE
% age_pd = fitlm(data_pd.ENROLL_AGE,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
% age_hc = fitlm(data_hc.ENROLL_AGE,NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC);
% 
% figure
% plot(age_pd)
% figure
% plot(age_hc)
% 
% 
% %% WEIGHT
% weight_pd = fitlm(data_pd.WGTKG,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
% weight_hc = fitlm(data_hc.WGTKG,NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC);
% 
% figure
% plot(weight_pd)
% figure
% plot(weight_hc)
% 
% %% MCATOT
% 
% mcatot_pd = fitlm(data_pd.MCATOT,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
% mcatot_hc = fitlm(data_hc.MCATOT,NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC);
% 
% figure
% plot(mcatot_hc)
% 
% figure
% plot(mcatot_pd)

% %% Prova lasso
% 
% figure(30)
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% subplot(131)
% coeff_model_caud_pd = lasso(regressors_pd_array.CAUDATE,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
% lasso_model_caud_pd = fitlm(regressors_pd_array.CAUDATE*coeff_model_caud_pd,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
% plot(lasso_model_caud_pd)
% % xlim([-0.5 0.6])
% % ylim([-0.6 0.8])
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% title('Caudate linear fit - PD')
% 
% subplot(132)
% coeff_model_put_pd = lasso(regressors_pd_array.PUTAMEN,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD);
% lasso_model_put_pd = fitlm(regressors_pd_array.PUTAMEN*coeff_model_put_pd,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD);
% 
% plot(lasso_model_put_pd)
% % xlim([-0.5 0.6])
% % ylim([-0.6 0.8])
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% title('Putamen linear fit - PD')
% 
% subplot(133)
% coeff_model_put_ant_pd = lasso(regressors_pd_array.PUTAMEN_ANT,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD);
% lasso_model_put_ant_pd = fitlm(regressors_pd_array.PUTAMEN_ANT*coeff_model_put_ant_pd,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD);
% 
% plot(lasso_model_put_ant_pd)
% % xlim([-0.5 0.6])
% % ylim([-0.6 0.8])
% xlabel('Covariates and symptoms')
% ylabel('Lateralization index')
% title('Putamen Anterior linear fit - PD')
% sgtitle('Lateralization index - Regressors FIT PD')










%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CESTINO
% % only NP1R, NP3, NP4 can have 101
% idx_tests = [35:40, 61:93, 98:103];
% for i = 1:length(idx_tests)
%     idx_samples.npTEST_101.(variables{idx_tests(i)}) = find(data_pd.(variables{idx_tests(i)}) == 101);
% end
% % Remove from struct empty fields
% fields = fieldnames(idx_samples.npTEST_101);
% idx_samples.npTEST_101 = rmfield(idx_samples.npTEST_101, fields(structfun(@isempty, idx_samples.npTEST_101)));
% % Remove 101

%while model_caud_pd.Rsquared.Ordinary < 0.2

% for i =1:size(covariates_pd,2)
%     for j =1:size(covariates_pd,2)
%          a = covariates_pd(:,i:j);
%          model_caud_pd = fitlm(covariates_pd,LATERALIZATION_coeff.CAUDATE.PD);
%          if model_caud_pd.Rsquared.Ordinary > 0.2
%              disp([i,' ',j,' ', model_caud_pd.Rsquared.Ordinary])
%          end
%     end
% end
% %end



