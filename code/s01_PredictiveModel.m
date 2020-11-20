% This scripts assumes you are running it from the code/ folder of the
% cloned github repository. Do not change directory structure otherwise you
% will need to update paths. 

%% Load saved variables
hcp = load('../data/HCP_icasig.mat');
abcd = load('../data/ABCD_icasig.mat');
load('../data/InSample.mat')

xt_hcp = bsxfun(@minus,featuremat,hcp.mu);
xt_abcd = bsxfun(@minus,featuremat,abcd.mu);

Abig_hcp = (pinv(hcp.icasig')*(xt_hcp)')';
Abig_abcd = (pinv(abcd.icasig')*(xt_abcd)')';

%% ABCD & HCP Top k components Predictive model

NumComp = 10;
NumComps = size(NumComp,2);
p = size(pheno,2);
pheno_pred_hcp = zeros(n,p,NumComps);
pheno_pred_abcd = zeros(n,p,NumComps);


    for iFold = 1:nFold
        
        %find train and test data for this fold
        test_idx = folds==iFold;
        train_idx = ~test_idx;
        n_train = sum(train_idx);
        n_test = sum(test_idx);
        k=NumComp;
        
        fprintf(1,'HCP & ABCD top k=10 PCA Fold %d of %d\n',iFold,nFold);

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
        
    end
    
% Get AUC ABCD & HCP top k
[x_pos_hcp,y_pos_hcp,t_pos_hcp,auc_pos_hcp] = perfcurve(pheno, pheno_pred_hcp,2);
[x_pos_abcd,y_pos_abcd,t_pos_abcd,auc_pos_abcd] = perfcurve(pheno, pheno_pred_abcd,2);
   
%% In sample Top k components Predictive model

pheno_pred_insamp = zeros(n,p,NumComps);

for iFold = 1:nFold
     fprintf(1,'.');
     %find train and test data for this fold
     test_idx = folds==iFold;
     train_idx = ~test_idx;
     n_train = sum(train_idx);
     n_test = sum(test_idx);
     k=NumComp;
     
    fprintf(1,'InSample top k=10 PCA Fold %d of %d\n',iFold,nFold);
    
    coeff = pca(featuremat(train_idx,:));
    components{iFold} = coeff';

    %mean center train, and mean center test with train means
    mu = mean(featuremat(train_idx,:));
    x = bsxfun(@minus,featuremat,mu);

    %calculate expressions for each subject for train and test
    A{iFold} = (pinv(components{iFold}')*x')';
    n_train = sum(train_idx);
    n_test = sum(test_idx);
    Abig = A{iFold}(train_idx,1:k);
    Abig_test = A{iFold}(test_idx,1:k);

    X = [ones(n_train,1) Abig nuisance(train_idx,:)];
    b = pinv(X'*X)*X'*pheno(train_idx,:);
    
    fit1 = mnrfit([Abig nuisance(train_idx,:)],pheno(train_idx,1));
    temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
    pheno_pred_insamp(test_idx,1) = temp(:,2);

end

[x_pos_insamp,y_pos_insamp,t_pos_insamp,auc_pos_insamp] = perfcurve(pheno, pheno_pred_insamp,2);

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
    pheno_pred_abcd2 = zeros(n,p,NumComps);
    for iFold = 1:nFold
        fprintf(1,'ABCD pheno-corr PCA Fold %d of %d\n',iFold,nFold);

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
        pheno_pred_abcd2(test_idx,1) = temp_abcd(:,2);
    end
fprintf(1,'\n');

[x_pos_abcd_pheno(:,iPheno),y_pos_abcd_pheno(:,iPheno),t_pos_abcd_pheno(:,iPheno),auc_pos_abcd_pheno(iPheno)] = perfcurve(pheno, pheno_pred_abcd2,2);
end

%% HCP pheno select components predictive model

for iPheno = 1:8
    NumComp = phenocomp_hcp{iPheno};
    NumComps = size(NumComp,2);
    pheno_pred_hcp2 = zeros(n,p,NumComps);
    for iFold = 1:nFold
        fprintf(1,'HCP pheno-corr PCA Fold %d of %d\n',iFold,nFold);
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
        pheno_pred_hcp2(test_idx,1) = temp_hcp(:,2);

    end
fprintf(1,'\n');

[x_pos_hcp_pheno(:,iPheno),y_pos_hcp_pheno(:,iPheno),t_pos_hcp_pheno(:,iPheno),auc_pos_hcp_pheno(iPheno)] = perfcurve(pheno, pheno_pred_hcp2,2);
end

%% Plot AUCROC curve for all models
figure; 
hold on
plot(x_pos_hcp, y_pos_hcp,'-r');
plot(x_pos_abcd, y_pos_abcd,'-b');
plot(x_pos_insamp,y_pos_insamp,'g-');
plot(x_pos_abcd_pheno(:,7), y_pos_abcd_pheno(:,7),'-m')
plot(x_pos_hcp_pheno(:,7), y_pos_hcp_pheno(:,7),'-c')
x_dummy = [0:1];
y_dummy = [0:1];
plot(x_dummy,y_dummy,':k')
hold off
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Whole Brain ROC Curve')
legend({'HCP top10','ABCD top10','InSamp top10','ABCD fluid intel','HCP fluid intel','Random chance'},'Location','southeast');

%% Fit linear model with all data (no CV, for visualization purposes)

% In-sample pca
coeffall = pca(featuremat);
icasig = coeffall';
mu = mean(featuremat);
x = bsxfun(@minus,featuremat,mu);
Aall = (pinv(icasig')*x')';

% HCP top10 consensus map model fit
NumComp = 10;
tmp1 = fitglm([Abig_hcp(:,1:NumComp) nuisance],pheno-1,'Distribution','binomial');
p1 = table2array(tmp1.Coefficients(2:end,4));
b1 = table2array(tmp1.Coefficients(2:end,1));
% ABCD top10 consensus map model fit
tmp2 = fitglm([Abig_abcd(:,1:NumComp) nuisance],pheno-1,'Distribution','binomial');
p2 = table2array(tmp2.Coefficients(2:end,4));
b2 = table2array(tmp2.Coefficients(2:end,1));
% In-sample top10 consensus map model fit
tmp3 = fitglm([Aall(:,1:NumComp) nuisance],pheno-1,'Distribution','binomial');
p3 = table2array(tmp3.Coefficients(2:end,4));
b3 = table2array(tmp3.Coefficients(2:end,1));
% HCP fluid intelligence consensus map model fit
NumComp = hcp.sig_idx_PMAT24_A_CR;
tmp4 = fitglm([Abig_hcp(:,NumComp) nuisance],pheno-1,'Distribution','binomial');
p4 = table2array(tmp4.Coefficients(2:end,4));
b4 = table2array(tmp4.Coefficients(2:end,1));
% ABCD fluid intelligence consensus map model fit
NumComp = abcd.sig_nihtbx_fluidcomp_uncorrected;
tmp5 = fitglm([Abig_abcd(:,NumComp) nuisance],pheno-1,'Distribution','binomial');
p5 = table2array(tmp5.Coefficients(2:end,4));
b5 = table2array(tmp5.Coefficients(2:end,1));

% Plot consensus maps 
cons10_hcp = zscore(sum(bsxfun(@times,hcp.icasig(1:10,:),b1(1:end-4)),1)')';
openfig('../figures/HCP_top10_consensus.fig','visible')

cons10_abcd = zscore(sum(bsxfun(@times,abcd.icasig(1:10,:),b2(1:end-4)),1)')';
openfig('../figures/ABCD_top10_consensus.fig','visible')

cons10_insamp = zscore(sum(bsxfun(@times,icasig(1:10,:),b3(1:end-4)),1)')';
openfig('../figures/InSample_top10_consensus.fig','visible')

cons_fluid_hcp = zscore(sum(bsxfun(@times,hcp.icasig(hcp.sig_idx_PMAT24_A_CR,:),b4(1:end-4)),1)')';
openfig('../figures/HCP_fluidintel_consensus.fig','visible')

cons_fluid_abcd = zscore(sum(bsxfun(@times,abcd.icasig(abcd.sig_nihtbx_fluidcomp_uncorrected,:),b5(1:end-4)),1)')';
openfig('../figures/ABCD_fluidintel_consensus.fig','visible')
