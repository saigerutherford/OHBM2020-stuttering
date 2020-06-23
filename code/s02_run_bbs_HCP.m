[coeff, score, latent, ~, exp] = pca(featuremat_hcp);
Abig = score;
icasig = coeff';

%mean center 
mu = mean(featuremat_hcp);
x = bsxfun(@minus,featuremat_hcp,mu);

%calculate expressions for each subject
A = (pinv(icasig')*x')';
    
CardSort_Unadj = zeros(909,2); %1
Flanker_Unadj = zeros(909,2); %2
PicVocab_Unadj = zeros(909,2); %3
ProcSpeed_Unadj = zeros(909,2); %4
PicSeq_Unadj = zeros(909,2); %5
Reading_Unadj = zeros(909,2); %6
ListSort_Unadj = zeros(909,2); %7
PMAT24_A_CR = zeros(909,2); %8
SCPT_TN = zeros(909,2); %9
SCPT_TP = zeros(909,2); %10
SCPT_SPEC = zeros(909,2); %11
SCPT_SEN = zeros(909,2); %12

for i = 1:909
    [CardSort_Unadj(i,1) CardSort_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,1));
    [Flanker_Unadj(i,1) Flanker_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,2));
    [PicVocab_Unadj(i,1) PicVocab_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,3));
    [ProcSpeed_Unadj(i,1) ProcSpeed_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,4));
    [PicSeq_Unadj(i,1) PicSeq_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,5));
    [Reading_Unadj(i,1) Reading_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,6));
    [ListSort_Unadj(i,1) ListSort_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,7));
    [PMAT24_A_CR(i,1) PMAT24_A_CR(i,2)] = corr(A(:,i),pheno_hcp(:,8));
    [SCPT_TN(i,1) SCPT_TN(i,2)] = corr(A(:,i),pheno_hcp(:,9));
    [SCPT_TP(i,1) SCPT_TP(i,2)] = corr(A(:,i),pheno_hcp(:,10));
    [SCPT_SPEC(i,1) SCPT_SPEC(i,2)] = corr(A(:,i),pheno_hcp(:,11));
    [SCPT_SEN(i,1) SCPT_SEN(i,2)] = corr(A(:,i),pheno_hcp(:,12));
end

sig_idx_CardSort_Unadj = find(CardSort_Unadj(:,2)<0.01);
sig_idx_Flanker_Unadj = find(Flanker_Unadj(:,2)<0.01);
sig_idx_PicVocab_Unadj = find(PicVocab_Unadj(:,2)<0.01);
sig_idx_ProcSpeed_Unadj = find(ProcSpeed_Unadj(:,2)<0.01);
sig_idx_PicSeq_Unadj = find(PicSeq_Unadj(:,2)<0.01);
sig_idx_Reading_Unadj = find(Reading_Unadj(:,2)<0.01);
sig_idx_ListSort_Unadj = find(ListSort_Unadj(:,2)<0.01);
sig_idx_PMAT24_A_CR = find(PMAT24_A_CR(:,2)<0.01);
sig_idx_SCPT_TN = find(SCPT_TN(:,2)<0.01);
sig_idx_SCPT_TP = find(SCPT_TP(:,2)<0.01);
sig_idx_SCPT_SPEC = find(SCPT_SPEC(:,2)<0.01);
sig_idx_SCPT_SEN = find(SCPT_SEN(:,2)<0.01);
