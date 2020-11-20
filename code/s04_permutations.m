nPerm = 10000;
NumComp = 10;
shuf_idx = zeros(n,nPerm);
shufpheno = zeros(size(pheno));

for iPerm = 1:nPerm
    shuf_idx(:,iPerm) = randperm(n);
    shufpheno(:,iPerm) = pheno(shuf_idx(:,iPerm));
end

auc_perm_hcp = zeros(nPerm,1);
auc_perm_abcd = zeros(nPerm,1);
pool = parpool(12);

parfor iPerm = 1:nPerm
    fprintf(1,'%d\n',iPerm);
    pheno_pred_hcp = zeros(size(pheno));
    pheno_pred_abcd = zeros(size(pheno));

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
    
    [x_pos_hcp,y_pos_hcp,t_pos_hcp,auc_posp_hcp] = perfcurve(shufpheno(:,iPerm), pheno_pred_hcp,2);
    auc_perm_hcp(iPerm) = auc_posp_hcp;
    [x_pos_abcd,y_pos_abcd,t_pos_abcd,auc_posp_abcd] = perfcurve(shufpheno(:,iPerm), pheno_pred_abcd,2);
    auc_perm_abcd(iPerm) = auc_posp_abcd;
end

hcp_top10_perm_pval = sum(auc_perm_hcp>=auc_pos_hcp)/nPerm

abcd_top10_perm_pval = sum(auc_perm_abcd>=auc_pos_abcd)/nPerm
