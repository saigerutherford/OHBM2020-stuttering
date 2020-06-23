%% Add needed paths for loading connectomes & phenotype data

addpath('/usr/local/MethodsCore/matlabScripts/')

ParamTemplate = '/net/pepper/youth_longitudinal/FirstLevel/SooEun/temp/NF102_1/Power_parameters.mat';
CorrTemplate = '/net/pepper/youth_longitudinal/FirstLevel/SooEun/temp/[Subject]/Power_corr.mat';

% Path to phenotype file
PhenotypeFile = '/net/pepper/youth_longitudinal/Scripts/SooEun/LogReg//sooeunscanfile_visit1included_updated3.csv'; % NvsS1

% Load phenotype file
dat = readtable(PhenotypeFile);
subs = dat.Subject;
ParamPathCheck = struct('Template',ParamTemplate,'mode','check');
ParamPath = mc_GenPath(ParamPathCheck);
param = load(ParamPath);
roiMNI = param.parameters.rois.mni.coordinates;

%% Load stuttering connectomes
Nodes = 264;
Edges = (Nodes*(Nodes-1))/2;
n = size(subs,1);
p = (Nodes*(Nodes-1))/2;

featuremat = zeros(n,p);

%load and z-transform connectivity matrices
for iSubject = 1:n
    fprintf(1,'Loading subject %d of %d\n',iSubject,n);
    path = strrep(CorrTemplate,'[Subject]',subs{iSubject});
    r = load(path,'rMatrix');
    r = r.rMatrix;
    z = mc_FisherZ(r);
    z = mc_flatten_upper_triangle(z);
    featuremat(iSubject,:) = z;
end

% Check if any connectomes are empty, change NaN's to 0's
nan = sum(isnan(featuremat),2)==size(featuremat,2);
bad = nan;
featuremat = featuremat(~bad,:);
dat = dat(~bad,:);
featuremat(isnan(featuremat))=0;
pheno = [dat.Group];
pheno = pheno(~bad,:);
p = size(pheno,2);
%% Load saved HCP & ABCD components, set necessary variables
hcp = load('/net/pepper/youth_longitudinal/Scripts/SooEun/LogReg/hcp_icasig_phenos.mat');
abcd = load('/net/pepper/youth_longitudinal/Scripts/SooEun/LogReg/abcd_icasig_phenos.mat');

xt_hcp = bsxfun(@minus,featuremat,hcp.mu);
xt_abcd = bsxfun(@minus,featuremat,abcd.mu);

Abig_hcp = (pinv(hcp.icasig')*(xt_hcp)')';
Abig_abcd = (pinv(abcd.icasig')*(xt_abcd)')';

%% Train/Test split for 10-fold CV
n = size(dat,1);
nFold = 10;
replace = 1;
if (nFold==n)
    replace = 0;
end
rng('default');
rng(12345);
folds = randsample(1:nFold,n,replace);

%% Control for nuisance variables within train/test, Method 1
nuisance = [dat.meanFD dat.meanFD.^2 dat.age_mri_months dat.sex];
refnuisancerep = [0 0 1 0]; % set all to zero, except age (set to mean)


for iFold = 1:nFold
    test_idx = folds==iFold;
    train_idx = ~test_idx;
    
    n_train = sum(train_idx);
    n_test = sum(test_idx);
    for iN = 1:size(nuisance,2)
        if(refnuisancerep(iN)==0)
            refnuisance(test_idx,iN) = 0;
        else
            refnuisance(test_idx,iN) = mean(nuisance(train_idx,iN));
        end
    end
end

%% Control for nuisance variables by subtracting them out, Method 2
% X = [ones(n,1) nuisance];
% b = pinv(X'*X)*X'*pheno;
% X0 = X;
% X0(:,1) = 0;
% pheno_res = pheno - X0*b;
% 
% b = pinv(X'*X)*X'*Abig_hcp;
% Abig_hcp_res = Abig_hcp - X0*b;
% 
% b = pinv(X'*X)*X'*Abig_abcd;
% Abig_abcd_res = Abig_abcd - X0*b;


%% ABCD & HCP Top k components Predictive model

NumComp = 30;
NumComps = size(NumComp,2);

pheno_pred_hcp = zeros(n,p,NumComps);
pheno_pred_abcd = zeros(n,p,NumComps);

x_pos_hcp = zeros(122,size(NumComp,1));
y_pos_hcp = zeros(122,size(NumComp,1));
t_pos_hcp = zeros(122,size(NumComp,1));
auc_pos_hcp = zeros(size(NumComp));
x_pos_abcd = zeros(122,size(NumComp,1));
y_pos_abcd = zeros(122,size(NumComp,1));
t_pos_abcd = zeros(122,size(NumComp,1));
auc_pos_abcd = zeros(size(NumComp));

for iComp = 1:NumComp
    for iFold = 1:nFold
        
        %find train and test data for this fold
        test_idx = folds==iFold;
        train_idx = ~test_idx;
        n_train = sum(train_idx);
        n_test = sum(test_idx);
        k=iComp;
        
        Abig2_hcp = Abig_hcp(train_idx,1:k);
        Abig_test_hcp = Abig_hcp(test_idx,1:k);
        X_hcp = [ones(n_train,1) Abig2_hcp];
        b_hcp = pinv(X_hcp'*X_hcp)*X_hcp'*pheno(train_idx,1);
        
        Abig2_abcd = Abig_abcd(train_idx,1:k);
        Abig_test_abcd = Abig_abcd(test_idx,1:k);
        X_abcd = [ones(n_train,1) Abig2_abcd];
        b_abcd = pinv(X_abcd'*X_abcd)*X_abcd'*pheno(train_idx,1);
        
        fit1_hcp = mnrfit([Abig2_hcp nuisance(train_idx,:)],pheno(train_idx,1));
        temp_hcp = mnrval(fit1_hcp,[Abig_test_hcp refnuisance(test_idx,:)]);
        pheno_pred_hcp(test_idx,1) = temp_hcp(:,2);
        
        fit1_abcd = mnrfit([Abig2_abcd nuisance(train_idx,:)],pheno(train_idx,1));
        temp_abcd = mnrval(fit1_abcd,[Abig_test_abcd refnuisance(test_idx,:)]);
        pheno_pred_abcd(test_idx,1) = temp_abcd(:,2);
        [~,~,~,auc_pos_fold_hcp(iFold,iComp)] = perfcurve(pheno(test_idx),pheno_pred_hcp(test_idx),2);
        [~,~,~,auc_pos_fold_abcd(iFold,iComp)] = perfcurve(pheno(test_idx),pheno_pred_abcd(test_idx),2);
        
        
    end
    
    % Get AUC ABCD & HCP top k
    [x_pos_hcp(:,iComp),y_pos_hcp(:,iComp),t_pos_hcp(:,iComp),auc_pos_hcp(iComp)] = perfcurve(pheno, pheno_pred_hcp,2);
    [x_pos_abcd(:,iComp),y_pos_abcd(:,iComp),t_pos_abcd(:,iComp),auc_pos_abcd(iComp)] = perfcurve(pheno, pheno_pred_abcd,2);
    
end
mean_auc_hcp = smean(auc_pos_fold_hcp,1);
std_auc_hcp = std(auc_pos_fold_hcp,[],1);
ts = tinv(0.975,nFold-1);
ci_auc_hcp = ts.*std_auc_hcp./sqrt(nFold);

% Plot NumComp (used) vs. AUC
x_ax = 1:30;
figure; plot(x_ax,auc_pos_hcp,'r-');
hold on
plot(x_ax,auc_pos_abcd,'b-');

%% In sample Top k components Predictive model
[coeff, score, latent, ~, exp] = pca(featuremat);
Abig = score;
icasig = coeff';

pheno_pred_insamp = zeros(n,p,NumComps);
x_pos_insamp = zeros(122,size(NumComp,1));
y_pos_insamp = zeros(122,size(NumComp,1));
t_pos_insamp = zeros(122,size(NumComp,1));
auc_pos_insamp = zeros(size(NumComp));

for iComp = 1:NumComp
    for iFold = 1:nFold
        fprintf(1,'.');
        %find train and test data for this fold
        test_idx = folds==iFold;
        train_idx = ~test_idx;
        n_train = sum(train_idx);
        n_test = sum(test_idx);
        k=iComp;
        
        Abig2 = Abig(train_idx,1:k);
        Abig_test = Abig(test_idx,1:k);
        X = [ones(n_train,1) Abig2];
        b = pinv(X'*X)*X'*pheno(train_idx,1);
        
        fit1 = mnrfit([Abig2 nuisance(train_idx,:)],pheno(train_idx,1));
        temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
        pheno_pred_insamp(test_idx,1) = temp(:,2);
        
    end
    [x_pos_insamp(:,iComp),y_pos_insamp(:,iComp),t_pos_insamp(:,iComp),auc_pos_insamp(iComp)] = perfcurve(pheno, pheno_pred_insamp,2);
end

plot(x_ax,auc_pos_insamp,'g-');
% hold off
%% Define HCP & ABCD pheno-related components 
phenocomp_hcp = {hcp.sig_idx_CardSort_Unadj hcp.sig_idx_ListSort_Unadj hcp.sig_idx_PicVocab_Unadj ...
    hcp.sig_idx_ProcSpeed_Unadj hcp.sig_idx_PicSeq_Unadj hcp.sig_idx_Reading_Unadj hcp.sig_idx_PMAT24_A_CR hcp.sig_idx_SCPT_TP};

phenocomp_abcd = {abcd.sig_nihtbx_cardsort_uncorrected abcd.sig_nihtbx_list_uncorrected abcd.sig_nihtbx_picvocab_uncorrected ...
    abcd.sig_nihtbx_pattern_uncorrected abcd.sig_nihtbx_picture_uncorrected abcd.sig_nihtbx_reading_uncorrected ...
    abcd.sig_nihtbx_fluidcomp_uncorrected abcd.sig_tfmri_sst_all_beh_crgo_nt};

%% ABCD pheno components predictive model

for iPheno = 1:8
    NumComp = phenocomp_abcd{iPheno};
    NumComps = size(NumComp,2);
    pheno_pred_abcd = zeros(n,p,NumComps);
for iFold = 1:nFold
    fprintf(1,'.');
    %find train and test data for this fold
    test_idx = folds==iFold;
    train_idx = ~test_idx;
    n_train = sum(train_idx);
    n_test = sum(test_idx);
    k=NumComp;
    
    Abig2_abcd = Abig_abcd(train_idx,k);
    Abig_test_abcd = Abig_abcd(test_idx,k);
    X_abcd = [ones(n_train,1) Abig2_abcd];
    b_abcd = pinv(X_abcd'*X_abcd)*X_abcd'*pheno(train_idx,1);
    
    fit1_abcd = mnrfit([Abig2_abcd nuisance(train_idx,:)],pheno(train_idx,1));
    temp_abcd = mnrval(fit1_abcd,[Abig_test_abcd refnuisance(test_idx,:)]);
    pheno_pred_abcd(test_idx,1) = temp_abcd(:,2);
end
fprintf(1,'\n');

[x_pos_abcd,y_pos_abcd,t_pos_abcd,auc_pos_abcd_pheno(iPheno)] = perfcurve(pheno, pheno_pred_abcd,2);
end
% OHBM figure, only fluid intelligence & sustained attention
x_ax_pheno = [13 13];
plot(x_ax_pheno, auc_pos_abcd_pheno(7:8),'bo')
% All Phenos
% x_ax_pheno = [8 9 13 13 15 18 13 13];
% figure; plot(x_ax_pheno, auc_pos_abcd_pheno,'bo')
% hold on
%% HCP pheno select components predictive model

for iPheno = 1:8
    NumComp = phenocomp_hcp{iPheno};
    NumComps = size(NumComp,2);
    pheno_pred_hcp = zeros(n,p,NumComps);
for iFold = 1:nFold
    fprintf(1,'.');
    %find train and test data for this fold
    test_idx = folds==iFold;
    train_idx = ~test_idx;
    n_train = sum(train_idx);
    n_test = sum(test_idx);
    k=NumComp;
    
    Abig2_hcp = Abig_hcp(train_idx,k);
    Abig_test_hcp = Abig_hcp(test_idx,k);
    X_hcp = [ones(n_train,1) Abig2_hcp];
    b_hcp = pinv(X_hcp'*X_hcp)*X_hcp'*pheno(train_idx,1);
    
    fit1_hcp = mnrfit([Abig2_hcp nuisance(train_idx,:)],pheno(train_idx,1));
    temp_hcp = mnrval(fit1_hcp,[Abig_test_hcp refnuisance(test_idx,:)]);
    pheno_pred_hcp(test_idx,1) = temp_hcp(:,2);
    
end
fprintf(1,'\n');
[x_pos_hcp,y_pos_hcp,t_pos_hcp,auc_pos_hcp_pheno(iPheno)] = perfcurve(pheno, pheno_pred_hcp,2);
end

% OHBM poster only fluid intelligence & sustained attn
x_ax_pheno = [17 17];
plot(x_ax_pheno, auc_pos_hcp_pheno(7:8),'ro')
% All 
% x_ax_pheno = [14 21 17 15 22 18 17 17];
% plot(x_ax_pheno, auc_pos_hcp_pheno,'ro')
hold off
%% Fit linear model with all data (no CV, for visualization purposes)

NumComp = 5;
tmp1 = fitglm([Abig_hcp(:,1:NumComp) nuisance],pheno-1,'Distribution','binomial');
p1 = table2array(tmp1.Coefficients(2:end,4));
b1 = table2array(tmp1.Coefficients(2:end,1));

tmp2 = fitglm([Abig_abcd(:,1:NumComp) nuisance],pheno-1,'Distribution','binomial');
p2 = table2array(tmp2.Coefficients(2:end,4));
b2 = table2array(tmp2.Coefficients(2:end,1));

tmp3 = fitglm([Abig(:,1:NumComp) nuisance],pheno-1,'Distribution','binomial');
p3 = table2array(tmp3.Coefficients(2:end,4));
b3 = table2array(tmp3.Coefficients(2:end,1));

NumComp = hcp.sig_idx_PMAT24_A_CR;
tmp4 = fitglm([Abig_hcp(:,NumComp) nuisance],pheno-1,'Distribution','binomial');
p4 = table2array(tmp4.Coefficients(2:end,4));
b4 = table2array(tmp4.Coefficients(2:end,1));

NumComp = abcd.sig_nihtbx_fluidcomp_uncorrected;
tmp5 = fitglm([Abig_abcd(:,NumComp) nuisance],pheno-1,'Distribution','binomial');
p5 = table2array(tmp5.Coefficients(2:end,4));
b5 = table2array(tmp5.Coefficients(2:end,1));

% %% Visualize components
param = load(ParamTemplate);
roiMNI = param.parameters.rois.mni.coordinates;
NetsFile = '/net/parasite/HCP/Scripts/slab/PCA/FundDiffRepo/Data/Power_Nets.csv';

netscsv = dataset('File',NetsFile,'Delimiter',',');
netsmni = [netscsv.X netscsv.Y netscsv.Z netscsv.Network];
[a,b] = ismember(roiMNI,netsmni(:,1:3),'rows');
nets = netsmni(b,4);

% Plot upper triangle expression of components 1-5 for ABCD, HCP, Insample

cons5_hcp = zscore(sum(bsxfun(@times,hcp.icasig(1:5,:),b1(1:end-4)),1)')';
plot_jica_component(cons5_hcp,1,1,2,nets,'',[1:13]);

% mdl1 = fitglm([Abig_hcp(:,1:5) nuisance],pheno-1,'Distribution','binomial','Link','logit');
% mdl1 = fitglm([Abig_hcp(:,1:5) ],pheno-1,'Distribution','binomial','Link','logit');
% %tmp = [ones(34716,1) hcp.icasig(1:5,:)']*mdl1.Coefficients.Estimate;%(2:end);
% tmp = hcp.icasig(1:5,:)'*mdl1.Coefficients.Estimate(2:6);
% tmp = tmp - mean(tmp);
% tmp2 = featuremat*tmp;
% tmp3 = 1./(1+exp(-tmp2));
% corr(mdl1.Fitted.Probability,tmp3)
% plot_jica_component(tmp',1,1,2,nets,'',[1:13]);
% tmp4 = 1./(1+exp(-tmp));
% tmp5 = featuremat*(tmp4-mean(tmp4));
% corr(mdl1.Fitted.Probability,tmp5)
% plot_jica_component(tmp4',1,1,2,nets,'',[1:13]);
% 
% % pmap5_hcp = 1/(1+exp(-tmp));
% % pmap5_hcp = pmap5_hcp - mean(pmap5_hcp);
% % plot_jica_component(pmap5_hcp,1,1,2,nets,'',[1:13]);


cons5_abcd = zscore(sum(bsxfun(@times,abcd.icasig(1:5,:),b2(1:end-4)),1)')';
plot_jica_component(cons5_abcd,1,1,2,nets,'',[1:13]);

cons5_insamp = zscore(sum(bsxfun(@times,icasig(1:5,:),b3(1:end-4)),1)')';
plot_jica_component(cons5_insamp,1,1,2,nets,'',[1:13]);

cons_fluid_hcp = zscore(sum(bsxfun(@times,hcp.icasig(hcp.sig_idx_PMAT24_A_CR,:),b4(1:end-4)),1)')';
plot_jica_component(cons_fluid_hcp,1,1,2,nets,'',[1:13]);

cons_fluid_abcd = zscore(sum(bsxfun(@times,abcd.icasig(abcd.sig_nihtbx_fluidcomp_uncorrected,:),b5(1:end-4)),1)')';
plot_jica_component(cons_fluid_abcd,1,1,2,nets,'',[1:13]);