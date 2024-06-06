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


%%   Decide what to remove
% Remove from idx_nan empty fields -------- PD
fields_nan_pd = fieldnames(idx_nan.PD);
idx_nan.PD = rmfield(idx_nan.PD, fields_nan_pd(structfun(@isempty, idx_nan.PD)));

% Remove from idx_nan empty fields -------- HC
fields_nan_hc = fieldnames(idx_nan.HC);
idx_nan.HC = rmfield(idx_nan.HC, fields_nan_hc(structfun(@isempty, idx_nan.HC)));

%%    - Remove all PD with NaNs in some np test
idx_fields_tests = [26:54, 60:93, 97:103];
idx_nan_all_tests = [];
for i = 1:length(idx_fields_tests)
    idx_nan_all_tests = [idx_nan_all_tests; idx_nan.PD.(fields_nan_pd{idx_fields_tests(i)})];
end
idx_nan_all_tests = unique(idx_nan_all_tests);
clear i idx_fields_tests
% Sono 857 i pazienti a cui manca almeno una misura di un test 

% togliamo np4
idx_fields_tests_no_np4 = [26:54, 60:93];
idx_nan_tests_no_np4 = [];
for i = 1:length(idx_fields_tests_no_np4)
    idx_nan_tests_no_np4 = [idx_nan_tests_no_np4; idx_nan.PD.(fields_nan_pd{idx_fields_tests_no_np4(i)})];
end
idx_nan_tests_no_np4 = unique(idx_nan_tests_no_np4);
clear i idx_fields_tests_no_np4
% Sono 220

% Remove patients missing at least one measure of test, excluding np4
data_pd(idx_nan_tests_no_np4,:) = [];

%%    - TRY Remove all HC with NaNs in some np test
idx_fields_tests = [26:54, 60:93, 97:103]; % same indeces
idx_nan_all_tests = [];
for i = 1:length(idx_fields_tests)
    idx_nan_all_tests = [idx_nan_all_tests; idx_nan.HC.(fields_nan_pd{idx_fields_tests(i)})];
end
idx_nan_all_tests = unique(idx_nan_all_tests);
clear i idx_fields_tests
% Sono 253 i pazienti a cui manca almeno una misura di un test 

% togliamo np4
idx_fields_tests_no_np4 = [26:54, 60:93];
idx_nan_tests_no_np4 = [];
for i = 1:length(idx_fields_tests_no_np4)
    idx_nan_tests_no_np4 = [idx_nan_tests_no_np4; idx_nan.HC.(fields_nan_pd{idx_fields_tests_no_np4(i)})];
end
idx_nan_tests_no_np4 = unique(idx_nan_tests_no_np4);
clear i idx_fields_tests_no_np4
% Sono 37

% Remove hc missing at least one measure of test, excluding np4
data_hc(idx_nan_tests_no_np4,:) = [];

%%    - Remove patients without HAND indication
data_pd(idx_nan.PD.HANDED, :) = [];
% HC don't have nan in hand

%%    - Remove patients with low QUALITY RATING
data_pd(find(data_pd.DATSCAN_QUALITY_RATING==3),:) = [];
data_hc(find(data_pd.DATSCAN_QUALITY_RATING==3),:) = [];

%%    - Remove patients with not COMPLETED DATSCAN
data_pd(find(data_pd.DATSCAN==0),:) = [];
data_hc(find(data_pd.DATSCAN==0),:) = [];


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
    y = [data_pd.(variables{cont_variables(i)}); data_hc.(variables{cont_variables(i)})]';
    LAT_groups = [ones(1,size(data_pd.(variables{cont_variables(i)}),1)), 2*ones(1,size(data_hc.(variables{cont_variables(i)}),1))];
    anova_test.(variables{cont_variables(i)}).p = anova1(y, LAT_groups);
    
    if anova_test.(variables{cont_variables(i)}).p > 0.05
        disp(strcat(string(variables{cont_variables(i)}), ': Accepted null hyp (same mean)'))
    else 
        disp(strcat(string(variables{cont_variables(i)}), ': Rejected null hyp (different mean)'))
    end

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
            disp(strcat(string(variables{discr_variables(2)}), ': LAT', region, ' Accepted null hyp (same mean)'))
        else 
            disp(strcat(string(variables{discr_variables(2)}), ': LAT', region, ' Rejected null hyp (different mean)'))
        end
        
        clear y groups_anova
    
    end
    clear i j region
end


% Difference in PUTAMEN and PUTAMEN ANT lateralization between male and
% female
figure
boxplot(LAT_groups.HC.PUTAMEN.Male, 'Colors','r')
hold on
boxplot(LAT_groups.HC.PUTAMEN.Female,'Colors','b')

figure
boxplot(LAT_groups.HC.PUTAMEN_ANT.Male, 'Colors','r')
hold on
boxplot(LAT_groups.HC.PUTAMEN_ANT.Female,'Colors','b')


<<<<<<< Updated upstream
%% LINEAR REGRESSION of variables of interest
=======



%% CORRELATION MATRIX of variables of interest
>>>>>>> Stashed changes
% fare scatter
%%   - HC
idx_covariates = [4,153,154,34, 55, 94, 158,27:54, 60:93];
covariates_hc = data_hc(:, idx_covariates); % age, weight, height,  np test + mcatot
covariates_hc.('Caudate lat coeff') = LATERALIZATION_coeff.CAUDATE.HC;
covariates_hc.('Putamen lat coeff') = LATERALIZATION_coeff.PUTAMEN.HC;
covariates_hc.('Putamen ant lat coeff') = LATERALIZATION_coeff.PUTAMEN_ANT.HC;

[p_corr_hc, R_corr_matrix_hc] = corrcoef(table2array(covariates_hc), 'Rows', 'pairwise');

figure
imagesc(R_corr_matrix_hc)
colormap parula
colorbar
xticks(1:width(covariates_hc))
xticklabels(covariates_hc.Properties.VariableNames)
yticks(1:width(covariates_hc))
yticklabels(covariates_hc.Properties.VariableNames)
title("Correlation HC age, weight, height,  np test + mcatot")

covariates_to_save_hc = [];
for i =1:size(R_corr_matrix_hc,1)
    for j = 1:size(R_corr_matrix_hc,2)
        if p_corr_hc(i,j) < 0.05 && R_corr_matrix_hc(i,j) > 0.5
                if (j < (size(R_corr_matrix_hc,1) -3)) && (i > (size(R_corr_matrix_hc,1) -3)) 
                    disp([covariates_hc.Properties.VariableNames{i}, ' correlated with ', covariates_hc.Properties.VariableNames{j}])
                    covariates_to_save_hc = [covariates_to_save_hc,  convertCharsToStrings(covariates_hc.Properties.VariableNames{j})];
                end
       end
    end
end

covariates_to_save_hc = unique(covariates_to_save_hc);
%%  - PD
covariates_pd = data_pd(:, idx_covariates); % age, weight, height,  np test + mcatot
covariates_pd.('Caudate lat coeff') = LATERALIZATION_coeff.CAUDATE.PD;
covariates_pd.('Putamen lat coeff') = LATERALIZATION_coeff.PUTAMEN.PD;
covariates_pd.('Putamen ant lat coeff') = LATERALIZATION_coeff.PUTAMEN_ANT.PD;

[p_corr_pd, R_corr_matrix_pd] = corrcoef(table2array(covariates_pd), 'Rows', 'pairwise');

figure
imagesc(R_corr_matrix_pd)
colormap parula
colorbar
xticks(1:width(covariates_pd))
xticklabels(covariates_pd.Properties.VariableNames)
yticks(1:width(covariates_pd))
yticklabels(covariates_pd.Properties.VariableNames)
title("Correlation HC age, weight, height,  np test + mcatot")

covariates_to_save_pd = [];
for i =1:size(R_corr_matrix_pd,1)
    for j = 1:size(R_corr_matrix_pd,2)
        if p_corr_pd(i,j) < 0.05 && R_corr_matrix_pd(i,j) > 0.5
                if (j < (size(R_corr_matrix_pd,1) -3)) && (i > (size(R_corr_matrix_pd,1) -3)) 
                    disp([covariates_pd.Properties.VariableNames{i}, ' correlated with ', covariates_pd.Properties.VariableNames{j}])
                    covariates_to_save_pd = [covariates_to_save_pd,  convertCharsToStrings(covariates_pd.Properties.VariableNames{j})];
                end
       end
    end
end

covariates_to_save_pd  = unique(covariates_to_save_pd);

%% LINEAR REGRESSION
%% - HC
covariates_hc = table2array(covariates_hc(:,contains(covariates_hc.Properties.VariableNames,covariates_to_save_hc)));

figure
subplot(131)
model_caud_hc = fitlm(covariates_hc,NOT_ABS.LATERALIZATION_coeff.CAUDATE.HC);
plot(model_caud_hc)

subplot(132)
model_put_hc = fitlm(covariates_hc,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.HC);
plot(model_put_hc)

subplot(133)
model_put_ant_hc = fitlm(covariates_hc,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.HC);
plot(model_put_ant_hc)


% Statistics
%caudate
% anova
% anova_caud = anova(model_caud_hc);
% % coef test
% coed_test_caud = coefTest(model_caud_hc);
% % [pdep_caud,x,y] = partialDependence(model_caud,{'x1','x3'});
% % figure
% % imagesc(x,y,pdep_caud)
% % putamen
% % anova
% anova_put = anova(model_put_hc);
% coef test
% coed_test_caud = coefTest(model_put);
% [pdep_put,x,y] = partialDependence(model_put,{'x1','x3'});
% figure
% imagesc(x,y,pdep_put)

%% - PD
covariates_pd = table2array(covariates_pd(:,contains(covariates_pd.Properties.VariableNames,covariates_to_save_pd)));


figure
subplot(131)
model_caud_pd = fitlm(covariates_pd,NOT_ABS.LATERALIZATION_coeff.CAUDATE.PD);
plot(model_caud_pd)

subplot(132)
model_put_pd = fitlm(covariates_pd,NOT_ABS.LATERALIZATION_coeff.PUTAMEN.PD);
plot(model_put_pd)

subplot(133)
model_put_ant_pd = fitlm(covariates_pd,NOT_ABS.LATERALIZATION_coeff.PUTAMEN_ANT.PD);
plot(model_put_ant_pd)


% Statistics
%caudate
% anova
% anova_caud = anova(model_caud_l);
% coef test
% coed_test_caud = coefTest(model_caud_l);
% [pdep_caud,x,y] = partialDependence(model_caud,{'x1','x3'});
% figure
% imagesc(x,y,pdep_caud)
% putamen
% anova
% anova_put = anova(model_put_l);
% coef test
% coed_test_caud = coefTest(model_put);
% [pdep_put,x,y] = partialDependence(model_put,{'x1','x3'});
% figure
% imagesc(x,y,pdep_put)



















































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
