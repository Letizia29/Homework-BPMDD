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

%%%%%%%
% Da gestire PDSTATE
% AV serve? Solo DATSCAN? Aggiungere AV nella conversione/analisi o
% toglierlo definitivamente dal dataset
%%%%%%%

%%   Extraction of NaNs indeces
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

%% IMPORTANT VARIABLES EXTRACTION

%%   Dominant hand division
idx_samples.HAND.PD.right = find(strcmp(data_pd.HANDED, 'Right'));
idx_samples.HAND.PD.left = find(strcmp(data_pd.HANDED, 'Left'));
idx_samples.HAND.PD.mixed = find(strcmp(data_pd.HANDED, 'Mixed'));

idx_samples.HAND.HC.right = find(strcmp(data_hc.HANDED, 'Right'));
idx_samples.HAND.HC.left = find(strcmp(data_hc.HANDED, 'Left'));
idx_samples.HAND.HC.mixed = find(strcmp(data_hc.HANDED, 'Mixed'));

%%   Neurophysiological testing and assessment: test in 4 parts
% normalizzazione: per poter confrontare NPTOT, occorre togliere prima
% tutti i 101, poi dividere la somma per il numero di test fatti. Quindi
% fare una media per ogni NP1P, NP1R, NP2, NP3, NP4, escludendo i 101.
% MA i 101 sono già tolti

np_test_normalized.np1ptot = data_pd.NP1PTOT ./ length(27:33);
np_test_normalized.np1rtot = data_pd.NP1RTOT ./ length(35:40);
np_test_normalized.np2ptot = data_pd.NP2PTOT ./ length(42:54);
np_test_normalized.np3tot = data_pd.NP3TOT ./ length(61:93);
np_test_normalized.np4tot = data_pd.NP4TOT ./ length(98:103);

% plot per curiosità
figure, hold on
plot(np_test_normalized.np1ptot)
plot(np_test_normalized.np1rtot)
plot(np_test_normalized.np2ptot)
plot(np_test_normalized.np3tot)
plot(np_test_normalized.np4tot)
axis tight
xlabel('patients')
title('normalized NPTEST')
ylabel('test score')
legend('1p', '1r', '2p', '3', '4')
ylim([-0.5 4])

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
    groups = [ones(1,size(data_pd.(variables{cont_variables(i)}),1)), 2*ones(1,size(data_hc.(variables{cont_variables(i)}),1))];
    anova_test.(variables{cont_variables(i)}).p = anova1(y, groups);
    
    if anova_test.(variables{cont_variables(i)}).p > 0.05
        disp(strcat(string(variables{cont_variables(i)}), ': Accepted null hyp (same mean)'))
    else 
        disp(strcat(string(variables{cont_variables(i)}), ': Rejected null hyp (different mean)'))
    end

end

clear i y groups

%% DISCRETE VARIABLES
discr_variables = [7, 8, 9, 11]; % ethnicity, sex, hand, familiarity

% Conversion to numerical levels
for i = 1:length(discr_variables)-1 % no need for familiarity
    legend_pd.(variables{discr_variables(i)}) = unique(data_pd.(variables{discr_variables(i)}));
    for j = 1:length(legend_pd.(variables{discr_variables(i)}))
        ind_temp = find(strcmp(data_pd.(variables{discr_variables(i)}), legend_pd.(variables{discr_variables(i)})(j)));
        for k = 1:length(ind_temp)
            data_pd.(variables{discr_variables(i)}){ind_temp(k)} = j;
        end
    end
    data_pd.(variables{discr_variables(i)}) = cellfun(@double, (data_pd.(variables{discr_variables(i)})));
end

% Conversion to numerical levels
for i = 1:length(discr_variables)-1 % no need for familiarity
    legend_hc.(variables{discr_variables(i)}) = unique(data_hc.(variables{discr_variables(i)}));
    for j = 1:length(legend_hc.(variables{discr_variables(i)}))
        ind_temp = find(strcmp(data_hc.(variables{discr_variables(i)}), legend_hc.(variables{discr_variables(i)})(j)));
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

%% MANAGING OUTLIERS (Box plots)







%% LATERALIZATION DATA
% Remove not completed DAT SCANS ( DATSCAN == 0 ) % è da fare prima?????
idx_samples.PD.DATSCAN_incomplete = find(data_pd.DATSCAN == 0); % 11
data_pd(idx_samples.PD.DATSCAN_incomplete,:) = [];

idx_samples.HC.DATSCAN_incomplete  = find(data_hc.DATSCAN == 0);
data_hc(idx_samples.HC.DATSCAN_incomplete ,:) = [];

% Extraction of lateralization index
% Right = |(right - left)/(right + left)| > 0.20
% Left = |(left - right)/(right + left)| > 0.20
% From: Ipsilateral deficits of dopaminergic neurotransmission in Parkinson s disease

% DATSCAN lateralization PD ------------------------------------------
LATERALIZATION_index.CAUDATE.PD = (data_pd.DATSCAN_CAUDATE_R - data_pd.DATSCAN_CAUDATE_L)./(data_pd.DATSCAN_CAUDATE_R + data_pd.DATSCAN_CAUDATE_L);
LATERALIZATION_index.PUTAMEN.PD = (data_pd.DATSCAN_PUTAMEN_R - data_pd.DATSCAN_PUTAMEN_L)./(data_pd.DATSCAN_PUTAMEN_R + data_pd.DATSCAN_PUTAMEN_L);
LATERALIZATION_index.PUTAMEN_ANT.PD = (data_pd.DATSCAN_PUTAMEN_R_ANT - data_pd.DATSCAN_PUTAMEN_L_ANT)./(data_pd.DATSCAN_PUTAMEN_R_ANT + data_pd.DATSCAN_PUTAMEN_L_ANT);

% DATSCAN lateralization HC ------------------------------------------
LATERALIZATION_index.CAUDATE.HC = (data_hc.DATSCAN_CAUDATE_R - data_hc.DATSCAN_CAUDATE_L)./(data_hc.DATSCAN_CAUDATE_R + data_hc.DATSCAN_CAUDATE_L);
LATERALIZATION_index.PUTAMEN.HC = (data_hc.DATSCAN_PUTAMEN_R - data_hc.DATSCAN_PUTAMEN_L)./(data_hc.DATSCAN_PUTAMEN_R + data_hc.DATSCAN_PUTAMEN_L);
LATERALIZATION_index.PUTAMEN_ANT.HC = (data_hc.DATSCAN_PUTAMEN_R_ANT - data_hc.DATSCAN_PUTAMEN_L_ANT)./(data_hc.DATSCAN_PUTAMEN_R_ANT + data_hc.DATSCAN_PUTAMEN_L_ANT);

% Grafici distribuzioni lateralizzazione?

% If the difference right-left is >0, then the DAT signal from the right
% side of the ROI is stronger than the left side.

% With significant difference: 20%

% save indeces of relevant lateralizations - PD
for i = 1:length(ROIs_labels)
    idx_samples.LATERALIZATION.PD.(ROIs_labels(i)).right = find(LATERALIZATION_index.(ROIs_labels(i)).PD > 0.2);
    idx_samples.LATERALIZATION.PD.(ROIs_labels(i)).left = find(LATERALIZATION_index.(ROIs_labels(i)).PD < -0.2);
    idx_samples.LATERALIZATION.PD.(ROIs_labels(i)).none = find(LATERALIZATION_index.(ROIs_labels(i)).PD < 0.2 & LATERALIZATION_index.(ROIs_labels(i)).PD > -0.2 );
end

% save indeces of relevant lateralizations - HC
for i = 1:length(ROIs_labels)
    idx_samples.LATERALIZATION.HC.(ROIs_labels(i)).right = find(LATERALIZATION_index.(ROIs_labels(i)).HC > 0.2);
    idx_samples.LATERALIZATION.HC.(ROIs_labels(i)).left = find(LATERALIZATION_index.(ROIs_labels(i)).HC < -0.2);
    idx_samples.LATERALIZATION.HC.(ROIs_labels(i)).none = find(LATERALIZATION_index.(ROIs_labels(i)).HC < 0.2 & LATERALIZATION_index.(ROIs_labels(i)).HC > -0.2 );
end


%% STATISTICAL ANALYSIS ------------- DA SISTEMARE

%% Correlation matrices



%% One-way Anova
% % Control gaussianity
% rois_names = fieldnames(DATSCAN);
% for i=2:length(rois_names)
%     rois_names_diag = fieldnames(DATSCAN.(rois_names{i,1}));
%     for j =1:length(rois_names_diag)
%         lil = lillietest(DATSCAN.(rois_names{i}).(rois_names_diag{j}));
%         if lil == 1
%             disp([rois_names{i}; rois_names_diag{j}; "not guassian"])
%         end
%     end
% end
% %%Absulte value
% % Caudate
% y = [abs(DATSCAN.CAUDATE_lat.PD)' abs(DATSCAN.CAUDATE_lat.HC)']';
% group1 = [repmat("PD_caudate",length(DATSCAN.CAUDATE_lat.PD),1);repmat("HC_caudate",length(DATSCAN.CAUDATE_lat.HC),1)];
% % [p,tbl,stats]  = anova1(y,group1);
% % 
% % % Putamen
% y = [abs(DATSCAN.PUTAMEN_lat.PD)' abs(DATSCAN.PUTAMEN_lat.HC)']';
% group2 = [repmat("PD_putamen",length(DATSCAN.PUTAMEN_lat.PD),1);repmat("HC_putamen",length(DATSCAN.PUTAMEN_lat.HC),1)];
% [p,tbl,stats]  = anova1(y,group2);
% % 
% % % Putamen ANT
% y = [abs(DATSCAN.PUTAMEN_ANT_lat.PD)' abs(DATSCAN.PUTAMEN_ANT_lat.HC)']';
% group3 = [repmat("PD_putamen_ant",length(DATSCAN.PUTAMEN_ANT_lat.PD),1);repmat("HC_putamen_ant",length(DATSCAN.PUTAMEN_ANT_lat.HC),1)];
% [p,tbl,stats]  = anova1(y,group3);
% % 
% y = [abs(DATSCAN.CAUDATE_lat.PD)' abs(DATSCAN.CAUDATE_lat.HC)' abs(DATSCAN.PUTAMEN_lat.PD)' abs(DATSCAN.PUTAMEN_lat.HC)' abs(DATSCAN.PUTAMEN_ANT_lat.PD)'  abs(DATSCAN.PUTAMEN_ANT_lat.HC)']';
% group4 = [group1; group2; group3];
% [p,tbl,stats]  = anova1(y,group4);
% % 
% 
% %multcompare(stats)
% %%Left
% %abs(DATSCAN.PUTAMEN_lat.PD)' abs(DATSCAN.PUTAMEN_lat.HC)' abs(DATSCAN.PUTAMEN_ANT_lat.PD)'  abs(DATSCAN.PUTAMEN_ANT_lat.HC)'
% % y = [DATSCAN.CAUDATE_lat.PD(lateralization.CAUDATE.PD.left.index)' DATSCAN.CAUDATE_lat.HC(lateralization.CAUDATE.HC.left.index)']';
% % group5 = [repmat("PD_caudate_LEFT",length(lateralization.CAUDATE.PD.left.index),1);repmat("HC_caudate_LEFT",length(lateralization.CAUDATE.HC.left.index),1)];
% % [p,tbl,stats]  = anova1(y,group5);
% 
%% Scatterplot
% j = [11,18,32,66];
% for j=j
%     figure
%     for i=1:3
%         subplot(1,3,i)
%         scatter(table2array(new_data_pd(:,i)),table2array(new_data_pd(:,j)))
%         hold on
%         scatter(table2array(new_data_hc(:,i)), table2array(new_data_hc(:,j)))
%         xlabel('Lateralization')
%         ylabel('Symptoms')
%         legend('PD','HC')
%         hold off
%         title(new_data_pd.Properties.VariableNames{j}, new_data_pd.Properties.VariableNames{i})
%     end
% end
% 
%% LINEAR REGRESSION of variables of interest
% %% right
% figure
% subplot(131)
% symptoms_r = [data_pd.NP3RIGN,data_pd.NP3RIGRU,data_pd.NP3RIGRL,data_pd.NP3PTRMR,data_pd.NP3KTRMR];
% model_caud_r = fitlm(symptoms_r,DATSCAN.CAUDATE_lat.PD);
% plot(model_caud_r)
% subplot(132)
% 
% symptoms_r = [data_pd.NP3RIGN,data_pd.NP3RIGRU,data_pd.NP3RIGRL,data_pd.NP3PTRMR,data_pd.NP3KTRMR];
% model_put_ant_r = fitlm(symptoms_r,DATSCAN.PUTAMEN_ANT_lat.PD);
% plot(model_put_ant_r)
% subplot(133)
% 
% symptoms_r = [data_pd.NP3RIGN,data_pd.NP3RIGRU,data_pd.NP3RIGRL,data_pd.NP3PTRMR,data_pd.NP3KTRMR];
% model_put_r = fitlm(symptoms_r,DATSCAN.PUTAMEN_lat.PD);
% plot(model_put_r)
% 
% % Statistics
% %caudate
% % anova
% anova_caud = anova(model_caud_r);
% % coef test
% coed_test_caud = coefTest(model_caud_r);
% % [pdep_caud,x,y] = partialDependence(model_caud,{'x1','x3'});
% % figure
% % imagesc(x,y,pdep_caud)
% % putamen
% % anova
% anova_put = anova(model_put_r);
% % coef test
% % coed_test_caud = coefTest(model_put);
% % [pdep_put,x,y] = partialDependence(model_put,{'x1','x3'});
% % figure
% % imagesc(x,y,pdep_put)
% %% left
% figure
% subplot(131)
% symptoms_l = [data_pd.NP3RIGN,data_pd.NP3RIGLU,data_pd.NP3RIGLL,data_pd.NP3PTRML,data_pd.NP3KTRML];
% model_caud_l = fitlm(symptoms_l,DATSCAN.CAUDATE_lat.PD);
% plot(model_caud_l)
% 
% subplot(132)
% model_put_ant_l = fitlm(symptoms_l,DATSCAN.PUTAMEN_ANT_lat.PD);
% plot(model_put_ant_l)
% 
% subplot(133)
% model_put_l = fitlm(symptoms_l,DATSCAN.PUTAMEN_lat.PD);
% plot(model_put_l)
% 
% % Statistics
% %caudate
% % anova
% anova_caud = anova(model_caud_l);
% % coef test
% coed_test_caud = coefTest(model_caud_l);
% % [pdep_caud,x,y] = partialDependence(model_caud,{'x1','x3'});
% % figure
% % imagesc(x,y,pdep_caud)
% % putamen
% % anova
% anova_put = anova(model_put_l);
% % coef test
% % coed_test_caud = coefTest(model_put);
% % [pdep_put,x,y] = partialDependence(model_put,{'x1','x3'});
% % figure
% % imagesc(x,y,pdep_put)



















































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
