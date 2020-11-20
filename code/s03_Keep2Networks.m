%% HCP, ABCD, In Sample Keep Two

Nets = max(nets);
for iNet = 1:(Nets-1)
    for jNet = (iNet+1):Nets
        fprintf(1,'Nets %d,%d\n',iNet,jNet);
        netmask = zeros(Nodes,Nodes);
        netmask(nets==iNet,nets==iNet) = 1;
        netmask(nets==iNet,nets==jNet) = 1;
        netmask(nets==jNet,nets==iNet) = 1;
        netmask(nets==jNet,nets==jNet) = 1;
        netmask = mc_flatten_upper_triangle(netmask)>0;
        NumComp = 10;
        NumComps = size(NumComp,2);
        n = size(pheno,1);
        p = size(pheno,2);

        % ABCD top 10
        featuremat_abcd2 = featuremat_abcd(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_abcd2);
        Abig = score;
        icasig = coeff';
        %mean center 
        mu = mean(featuremat_abcd2);
        x = bsxfun(@minus,featuremat_abcd2,mu);
        pheno_pred_abcd = zeros(n,p,NumComps);
        for iFold = 1:nFold
            %find train and test data for this fold
            test_idx = folds==iFold;
            train_idx = ~test_idx;
            n_train = sum(train_idx);
            n_test = sum(test_idx);
            k=NumComp;
            Abig2 = Abig(train_idx,1:k);
            Abig_test = Abig(test_idx,1:k);
            X = [ones(n_train,1) Abig2];
            b = pinv(X'*X)*X'*pheno(train_idx,1); 
            fit1 = mnrfit([Abig2 nuisance(train_idx,:)],pheno(train_idx,1));
            temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
            pheno_pred_abcd(test_idx,1) = temp(:,2);
        end
        [x_pos_abcd,y_pos_abcd,t_pos_abcd,auc_pos_abcd_network10(iNet,jNet)] = perfcurve(pheno, pheno_pred_abcd,2);

        
        % HCP top 10
        featuremat_hcp2 = featuremat_hcp(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_hcp2);
        Abig = score;
        icasig = coeff';
        NumComp = 10;
        NumComps = size(NumComp,2);
        %mean center 
        mu = mean(featuremat_hcp2);
        x = bsxfun(@minus,featuremat_hcp2,mu);
        pheno_pred_hcp = zeros(n,p,NumComps);
        for iFold = 1:nFold
            %find train and test data for this fold
            test_idx = folds==iFold;
            train_idx = ~test_idx;
            n_train = sum(train_idx);
            n_test = sum(test_idx);
            k=NumComp;
            Abig2 = Abig(train_idx,1:k);
            Abig_test = Abig(test_idx,1:k);
            X = [ones(n_train,1) Abig2];
            b = pinv(X'*X)*X'*pheno(train_idx,1); 
            fit1 = mnrfit([Abig2 nuisance(train_idx,:)],pheno(train_idx,1));
            temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
            pheno_pred_hcp(test_idx,1) = temp(:,2);
        end
        [x_pos_hcp,y_pos_hcp,t_pos_hcp,auc_pos_hcp_network10(iNet,jNet)] = perfcurve(pheno, pheno_pred_hcp,2);
        
        
        %In Sample top 10
        featuremat_2 = featuremat(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_2);
        Abig = score;
        icasig = coeff';
        pheno_pred = zeros(n,p,NumComps);
        for iFold = 1:nFold
            %find train and test data for this fold
            test_idx = folds==iFold;
            train_idx = ~test_idx;
            n_train = sum(train_idx);
            n_test = sum(test_idx);
            k=NumComp;
            Abig2 = Abig(train_idx,1:k);
            Abig_test = Abig(test_idx,1:k);
            X = [ones(n_train,1) Abig2];
            b = pinv(X'*X)*X'*pheno(train_idx,1); 
            fit1 = mnrfit([Abig2 nuisance(train_idx,:)],pheno(train_idx,1));
            temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
            pheno_pred(test_idx,1) = temp(:,2);
        end
        [x_pos,y_pos,t_pos,auc_pos_network10(iNet,jNet)] = perfcurve(pheno, pheno_pred,2);
    end
    
end

%% ABCD & HCP fluid intelligence keep 2
for iNet = 1:(Nets-1)
    for jNet = (iNet+1):Nets
        fprintf(1,'Nets %d,%d\n',iNet,jNet);
        netmask = zeros(Nodes,Nodes);
        netmask(nets==iNet,nets==iNet) = 1;
        netmask(nets==iNet,nets==jNet) = 1;
        netmask(nets==jNet,nets==iNet) = 1;
        netmask(nets==jNet,nets==jNet) = 1;
        netmask = mc_flatten_upper_triangle(netmask)>0;
        
        
        %ABCD fluid intel
        featuremat_abcd2 = featuremat_abcd(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_abcd2);
        Abig = score;
        icasig = coeff';
        %mean center 
        mu = mean(featuremat_abcd2);
        x = bsxfun(@minus,featuremat_abcd2,mu);
        %calculate expressions for each subject
        A = (pinv(icasig')*x')';
        comp = size(A,2);
        nihtbx_fluidcomp_uncorrected = zeros(size(A,2),2);
            for i = 1:comp
                [nihtbx_fluidcomp_uncorrected(i,1) nihtbx_fluidcomp_uncorrected(i,2)] = corr(A(:,i),pheno_abcd(:,8));
            end
        sig_nihtbx_fluidcomp_uncorrected = find(nihtbx_fluidcomp_uncorrected(:,2)<0.01);
        NumComp = sig_nihtbx_fluidcomp_uncorrected;
        NumComps = size(NumComp,2);
        pheno_pred_abcd = zeros(n,p,NumComps);
        for iFold = 1:nFold
            %find train and test data for this fold
            test_idx = folds==iFold;
            train_idx = ~test_idx;
            n_train = sum(train_idx);
            n_test = sum(test_idx);
            k=NumComp;
            Abig2 = Abig(train_idx,1:k);
            Abig_test = Abig(test_idx,1:k);
            X = [ones(n_train,1) Abig2];
            b = pinv(X'*X)*X'*pheno(train_idx,1); 
            fit1 = mnrfit([Abig2 nuisance(train_idx,:)],pheno(train_idx,1));
            temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
            pheno_pred_abcd(test_idx,1) = temp(:,2);
        end
        [x_pos_abcd,y_pos_abcd,t_pos_abcd,auc_pos_abcd_fluid(iNet,jNet)] = perfcurve(pheno, pheno_pred_abcd,2);
 
        
        % HCP fluid intel
        featuremat_hcp2 = featuremat_hcp(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_hcp2);
        Abig = score;
        icasig = coeff';
        %mean center 
        mu = mean(featuremat_hcp2);
        x = bsxfun(@minus,featuremat_hcp2,mu);
        %calculate expressions for each subject
        A = (pinv(icasig')*x')';
        comp = size(A,2);
        PMAT24_A_CR = zeros(size(A,2),2);
            for i = 1:comp
                [PMAT24_A_CR(i,1) PMAT24_A_CR(i,2)] = corr(A(:,i),pheno_hcp(:,8));
            end
        sig_idx_PMAT24_A_CR = find(PMAT24_A_CR(:,2)<0.01);
        NumComp = sig_idx_PMAT24_A_CR;
        NumComps = size(NumComp,2);
        pheno_pred_hcp = zeros(n,p,NumComps);
        for iFold = 1:nFold
            %find train and test data for this fold
            test_idx = folds==iFold;
            train_idx = ~test_idx;
            n_train = sum(train_idx);
            n_test = sum(test_idx);
            k=NumComp;
            Abig2 = Abig(train_idx,1:k);
            Abig_test = Abig(test_idx,1:k);
            X = [ones(n_train,1) Abig2];
            b = pinv(X'*X)*X'*pheno(train_idx,1); 
            fit1 = mnrfit([Abig2 nuisance(train_idx,:)],pheno(train_idx,1));
            temp = mnrval(fit1,[Abig_test refnuisance(test_idx,:)]);
            pheno_pred_hcp(test_idx,1) = temp(:,2);
        end
        [x_pos_hcp,y_pos_hcp,t_pos_hcp,auc_pos_hcp_fluid(iNet,jNet)] = perfcurve(pheno, pheno_pred_hcp,2);
    end
end