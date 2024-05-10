clear all
close all
clc


descr = readstruct("MasterPatient_Legend.json","FileType","json");
data = readtable("Patient_Master.csv");

DAT_SCAN_SPECT = [data.DATSCAN_CAUDATE_R data.DATSCAN_CAUDATE_L data.DATSCAN_PUTAMEN_R  data.DATSCAN_PUTAMEN_L data.DATSCAN_PUTAMEN_R_ANT data.DATSCAN_PUTAMEN_L_ANT];
%DAT_SCAN_PET =  [data.AV133_RCAUD_S data.AV133_LCAUD_S data.AV133_RPUTANT_S data.AV133_RPUTPOST_S data.AV133_LPUTANT_S data.AV133_LPUTPOST_S];

idx_HC = find(string(data.COHORT)=='HC');
HC_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_HC,:);

idx_PD = find(string(data.COHORT)=='PD');
PD_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_PD,:);

idx_SWEDD = find(string(data.COHORT)=='SWEDD');
SWEDD_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_SWEDD,:);

idx_Prodromal = find(string(data.COHORT)=='Prodromal');
Prodromal_DAT_SCAN_SPECT = DAT_SCAN_SPECT(idx_Prodromal,:);

NON_HC_DAT_SCAN_SPECT = DAT_SCAN_SPECT([idx_Prodromal; idx_PD; idx_SWEDD],:);

for i=1:6
    subplot(2,3,i)
    hist(NON_HC_DAT_SCAN_SPECT(:,i))
    hold on 
    hist(HC_DAT_SCAN_SPECT(:,i))
    hold off
end

