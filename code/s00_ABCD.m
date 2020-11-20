%% Part 1 ABCD: Load data
AnalysisDir = '/net/pepper/ABCD/Scripts/slab/PCA/repo';
addpath([AnalysisDir '/Utils/']);

DataDir = '/net/pepper/ABCD/FirstLevel/';

PhenotypeFile = '/net/pepper/youth_longitudinal/Scripts/SooEun/LogReg/ABCD_rest_data_excludephilips_phenos_dropna.csv';
Nodes = 264;
Edges = (Nodes*(Nodes-1))/2;

OutputPath = [AnalysisDir '/Results/'];

Model = 's6_Power264_p50f0b_nonaggr_p35mask';
ParamTemplate = [DataDir 'sub-NDARINV00NPMHND/avg/' Model '/' Model '_parameters.mat'];
CorrTemplate = [DataDir '[Subject]/avg/' Model '/' Model '_corr.mat'];

dat = readtable(PhenotypeFile);
subs_abcd = dat.Subject;
pheno_abcd = [dat.nihtbx_cardsort_uncorrected dat.nihtbx_list_uncorrected dat.nihtbx_picvocab_uncorrected  dat.nihtbx_pattern_uncorrected dat.nihtbx_picture_uncorrected dat.nihtbx_reading_uncorrected dat.nihtbx_fluidcomp_uncorrected dat.tfmri_sst_all_beh_crgo_nt];

n = size(subs_abcd,1);
p = (Nodes*(Nodes-1))/2;

featuremat_abcd = zeros(n,p);

%load and z-transform connectivity matrices
for iSubject = 1:n
    fprintf(1,'Loading subject %d of %d\n',iSubject,n);
    path = strrep(CorrTemplate,'[Subject]',subs_abcd{iSubject});
    r = load(path,'rMatrix');
    r = r.rMatrix;
    z = mc_FisherZ(r);
    z = mc_flatten_upper_triangle(z);
    featuremat_abcd(iSubject,:) = z;
end

%check for empty subjects
bad = sum(isnan(featuremat_abcd),2)==size(featuremat_abcd,2);
sum(bad)
featuremat_abcd(bad,:) = [];
featuremat_abcd(isnan(featuremat_abcd)) = 0;

subs = subs(~bad);
pheno_abcd = pheno_abcd(~bad,:);

%% Part 2 ABCD: PCA, choose components significantly correlated with phenotypes, save 
[coeff, score, latent, ~, exp] = pca(featuremat_abcd);
Abig = score;
icasig = coeff';

%mean center 
mu = mean(featuremat_abcd);
x = bsxfun(@minus,featuremat_abcd,mu);

%calculate expressions for each subject
A = (pinv(icasig')*x')';
    
nihtbx_cardsort_uncorrected = zeros(1508,2);
nihtbx_list_uncorrected = zeros(1508,2);
nihtbx_picvocab_uncorrected = zeros(1508,2);
nihtbx_pattern_uncorrected = zeros(1508,2);
nihtbx_picture_uncorrected = zeros(1508,2);
nihtbx_reading_uncorrected = zeros(1508,2);
nihtbx_fluidcomp_uncorrected = zeros(1508,2);
tfmri_sst_all_beh_crgo_nt = zeros(1508,2);

for i = 1:1508
    [nihtbx_cardsort_uncorrected(i,1) nihtbx_cardsort_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,1));
    [nihtbx_list_uncorrected(i,1) nihtbx_list_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,2));
    [nihtbx_picvocab_uncorrected(i,1) nihtbx_picvocab_uncorrected(i,2) ] = corr(A(:,i),pheno_abcd(:,3));
    [nihtbx_pattern_uncorrected(i,1) nihtbx_pattern_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,4));
    [nihtbx_picture_uncorrected(i,1) nihtbx_picture_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,5));
    [nihtbx_reading_uncorrected(i,1) nihtbx_reading_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,6));
    [nihtbx_fluidcomp_uncorrected(i,1)  nihtbx_fluidcomp_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,7));
    [tfmri_sst_all_beh_crgo_nt(i,1) tfmri_sst_all_beh_crgo_nt(i,2)] = corr(A(:,i),pheno_abcd(:,8));
end

sig_nihtbx_cardsort_uncorrected = find(nihtbx_cardsort_uncorrected(:,2)<0.01);
sig_nihtbx_list_uncorrected = find(nihtbx_list_uncorrected(:,2)<0.01);
sig_nihtbx_picvocab_uncorrected = find(nihtbx_picvocab_uncorrected(:,2)<0.01);
sig_nihtbx_pattern_uncorrected = find(nihtbx_pattern_uncorrected(:,2)<0.01);
sig_nihtbx_picture_uncorrected = find(nihtbx_picture_uncorrected(:,2)<0.01);
sig_nihtbx_reading_uncorrected = find(nihtbx_reading_uncorrected(:,2)<0.01);
sig_nihtbx_fluidcomp_uncorrected = find(nihtbx_fluidcomp_uncorrected(:,2)<0.01);
sig_tfmri_sst_all_beh_crgo_nt = find(tfmri_sst_all_beh_crgo_nt(:,2)<0.01);

% Clear unused variables & save workspace to be loaded in
% HCP_ABCD_CWS_RestingState_LogReg.m script

clear A AnalysisDir ans bad coeff CorrTemplate dat DataDir Edges exp featuremat_abcd ...
    i iSubject latent Model n nihtbx_cardsort_uncorrected nihtbx_fluidcomp_uncorrected ...
    nihtbx_list_uncorrected nihtbx_pattern_uncorrected nihtbx_picture_uncorrected ...
    nihtbx_picvocab_uncorrected nihtbx_reading_uncorrected Nodes OutputPath ...
    p ParamTemplate path pheno_abcd PhenotypeFile r score subs subs_abcd ...
    tfmri_sst_all_beh_crgo_nt x z

save('ABCD_icasig.mat')