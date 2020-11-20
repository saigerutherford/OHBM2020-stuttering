%% Part 1 HCP: Load data
addpath('/net/parasite/HCP/Scripts/slab/PCA/FundDiffRepo/Utils/');
DataDir = '/net/parasite/HCP/derivatives/FundDiffRepo/';
PhenotypeFile = '/net/pepper/youth_longitudinal/Scripts/SooEun/LogReg/HCP_unrestricted_phenotypes_include_addphenos.csv';

Nodes = 264;
Edges = (Nodes*(Nodes-1))/2;

OutputPath = '/net/parasite/HCP/Scripts/slab/PCA/FundDiffRepo/Results/';

Model = 's6_Power264_p50f0b_nonaggr_p35mask';
ParamTemplate = [DataDir 'Connectomes/100206/power_264/parameters.mat'];
CorrTemplate = [DataDir 'Connectomes/[Subject]/power_264/corr.mat'];

dat = readtable(PhenotypeFile);
subs_hcp = num2str(dat.Subject);

pheno_hcp = [dat.CardSort_Unadj	dat.ListSort_Unadj dat.PicVocab_Unadj dat.ProcSpeed_Unadj dat.PicSeq_Unadj dat.ReadEng_Unadj dat.PMAT24_A_CR dat.SCPT_TP];

n = size(subs_hcp,1);
p = (Nodes*(Nodes-1))/2;

featuremat_hcp = load_connectomes(CorrTemplate,subs_hcp);

%% Part 2 HCP: PCA, choose components significantly correlated with phenotypes, save 
[coeff, score, latent, ~, exp] = pca(featuremat_hcp);
Abig = score;
icasig = coeff';

%mean center 
mu = mean(featuremat_hcp);
x = bsxfun(@minus,featuremat_hcp,mu);

%calculate expressions for each subject
A = (pinv(icasig')*x')';
    
CardSort_Unadj = zeros(909,2);
ListSort_Unadj = zeros(909,2);
PicVocab_Unadj = zeros(909,2);
ProcSpeed_Unadj = zeros(909,2);
PicSeq_Unadj = zeros(909,2);
Reading_Unadj = zeros(909,2);
PMAT24_A_CR = zeros(909,2);
SCPT_TP = zeros(909,2);

for i = 1:909
    [CardSort_Unadj(i,1) CardSort_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,1));
    [ListSort_Unadj(i,1) ListSort_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,2));
    [PicVocab_Unadj(i,1) PicVocab_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,3));
    [ProcSpeed_Unadj(i,1) ProcSpeed_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,4));
    [PicSeq_Unadj(i,1) PicSeq_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,5));
    [Reading_Unadj(i,1) Reading_Unadj(i,2)] = corr(A(:,i),pheno_hcp(:,6));
    [PMAT24_A_CR(i,1) PMAT24_A_CR(i,2)] = corr(A(:,i),pheno_hcp(:,7));
    [SCPT_TP(i,1) SCPT_TP(i,2)] = corr(A(:,i),pheno_hcp(:,8));
end

sig_idx_CardSort_Unadj = find(CardSort_Unadj(:,2)<0.01);
sig_idx_ListSort_Unadj = find(ListSort_Unadj(:,2)<0.01);
sig_idx_PicVocab_Unadj = find(PicVocab_Unadj(:,2)<0.01);
sig_idx_ProcSpeed_Unadj = find(ProcSpeed_Unadj(:,2)<0.01);
sig_idx_PicSeq_Unadj = find(PicSeq_Unadj(:,2)<0.01);
sig_idx_Reading_Unadj = find(Reading_Unadj(:,2)<0.01);
sig_idx_PMAT24_A_CR = find(PMAT24_A_CR(:,2)<0.01);
sig_idx_SCPT_TP = find(SCPT_TP(:,2)<0.01);

% Clear unused variable and save workspace
clear A Abig CardSort_Unadj coeff CorrTemplate dat DataDir Edges exp featuremat_hcp ...
    i latent ListSort_Unadj Model n Nodes OutputPath p ParamTemplate pheno_hcp ...
    PhenotypeFile PicSeq_Unadj PicVocab_Unadj PMAT24_A_CR ProcSpeed_Unadj ... 
    Reading_Unadj score SCPT_TP subs_hcp x

save('HCP_icasig.mat')