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
pheno_abcd = [dat.BISBAS_fac dat.nihtbx_picvocab_uncorrected dat.nihtbx_flanker_uncorrected dat.nihtbx_list_uncorrected dat.nihtbx_cardsort_uncorrected dat.nihtbx_pattern_uncorrected dat.nihtbx_picture_uncorrected dat.nihtbx_reading_uncorrected dat.nihtbx_fluidcomp_uncorrected dat.nihtbx_cryst_uncorrected dat.tfmri_sst_all_beh_crgo_nt dat.tfmri_sst_all_beh_crlg_nt dat.tfmri_sst_all_beh_nrgo_nt dat.tfmri_sst_all_beh_crs_nt];

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
