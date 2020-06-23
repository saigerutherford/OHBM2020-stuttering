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

pheno_hcp = [dat.CardSort_Unadj	dat.Flanker_Unadj dat.PicVocab_Unadj dat.ProcSpeed_Unadj dat.PicSeq_Unadj dat.ReadEng_Unadj dat.ListSort_Unadj dat.PMAT24_A_CR dat.SCPT_TN dat.SCPT_TP dat.SCPT_SPEC dat.SCPT_SEN];

n = size(subs_hcp,1);
p = (Nodes*(Nodes-1))/2;

featuremat_hcp = load_connectomes(CorrTemplate,subs_hcp);