
[coeff, score, latent, ~, exp] = pca(featuremat_abcd);
Abig = score;
icasig = coeff';

%mean center 
mu = mean(featuremat_abcd);
x = bsxfun(@minus,featuremat_abcd,mu);

%calculate expressions for each subject
A = (pinv(icasig')*x')';
    

BISBAS_fac = zeros(1508,2);
nihtbx_picvocab_uncorrected = zeros(1508,2);
nihtbx_flanker_uncorrected = zeros(1508,2);
nihtbx_list_uncorrected = zeros(1508,2);
nihtbx_cardsort_uncorrected = zeros(1508,2);
nihtbx_pattern_uncorrected = zeros(1508,2);
nihtbx_picture_uncorrected = zeros(1508,2);
nihtbx_reading_uncorrected = zeros(1508,2);
nihtbx_fluidcomp_uncorrected = zeros(1508,2);
nihtbx_cryst_uncorrected = zeros(1508,2);
tfmri_sst_all_beh_crgo_nt = zeros(1508,2);
tfmri_sst_all_beh_crlg_nt = zeros(1508,2);
tfmri_sst_all_beh_nrgo_nt = zeros(1508,2);
tfmri_sst_all_beh_crs_nt = zeros(1508,2);


for i = 1:1508
    [BISBAS_fac(i,1) BISBAS_fac(i,2)] = corr(A(:,i),pheno_abcd(:,1));
    [nihtbx_picvocab_uncorrected(i,1) nihtbx_picvocab_uncorrected(i,2) ] = corr(A(:,i),pheno_abcd(:,2));
    [nihtbx_flanker_uncorrected(i,1) nihtbx_flanker_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,3));
    [nihtbx_list_uncorrected(i,1) nihtbx_list_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,4));
    [nihtbx_cardsort_uncorrected(i,1) nihtbx_cardsort_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,5));
    [nihtbx_pattern_uncorrected(i,1) nihtbx_pattern_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,6));
    [nihtbx_picture_uncorrected(i,1) nihtbx_picture_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,7));
    [nihtbx_reading_uncorrected(i,1) nihtbx_reading_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,8));
    [nihtbx_fluidcomp_uncorrected(i,1)  nihtbx_fluidcomp_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,9));
    [nihtbx_cryst_uncorrected(i,1) nihtbx_cryst_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,10));
    [tfmri_sst_all_beh_crgo_nt(i,1) tfmri_sst_all_beh_crgo_nt(i,2)] = corr(A(:,i),pheno_abcd(:,11));
    [tfmri_sst_all_beh_crlg_nt(i,1) tfmri_sst_all_beh_crlg_nt(i,2)] = corr(A(:,i),pheno_abcd(:,12));
    [tfmri_sst_all_beh_nrgo_nt(i,1) tfmri_sst_all_beh_nrgo_nt(i,2)] = corr(A(:,i),pheno_abcd(:,13));
    [tfmri_sst_all_beh_crs_nt(i,1) tfmri_sst_all_beh_crs_nt(i,2)] = corr(A(:,i),pheno_abcd(:,14));
end


sig_BISBAS_fac = find(BISBAS_fac(:,2)<0.01);
sig_nihtbx_picvocab_uncorrected = find(nihtbx_picvocab_uncorrected(:,2)<0.01);
sig_nihtbx_flanker_uncorrected = find(nihtbx_flanker_uncorrected(:,2)<0.01);
sig_nihtbx_list_uncorrected = find(nihtbx_list_uncorrected(:,2)<0.01);
sig_nihtbx_cardsort_uncorrected = find(nihtbx_cardsort_uncorrected(:,2)<0.01);
sig_nihtbx_pattern_uncorrected = find(nihtbx_pattern_uncorrected(:,2)<0.01);
sig_nihtbx_picture_uncorrected = find(nihtbx_picture_uncorrected(:,2)<0.01);
sig_nihtbx_reading_uncorrected = find(nihtbx_reading_uncorrected(:,2)<0.01);
sig_nihtbx_fluidcomp_uncorrected = find(nihtbx_fluidcomp_uncorrected(:,2)<0.01);
sig_nihtbx_cryst_uncorrected = find(nihtbx_cryst_uncorrected(:,2)<0.01);
sig_tfmri_sst_all_beh_crgo_nt = find(tfmri_sst_all_beh_crgo_nt(:,2)<0.01);
sig_tfmri_sst_all_beh_crlg_nt = find(tfmri_sst_all_beh_crlg_nt(:,2)<0.01);
sig_tfmri_sst_all_beh_nrgo_nt = find(tfmri_sst_all_beh_nrgo_nt(:,2)<0.01);
sig_tfmri_sst_all_beh_crs_nt = find(tfmri_sst_all_beh_crs_nt(:,2)<0.01);
