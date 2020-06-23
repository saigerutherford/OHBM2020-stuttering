addpath('/usr/local/MethodsCore/matlabScripts/')

clear auc_pos_abcd auc_pos_hcp auc_pos

ParamTemplate = '/net/pepper/youth_longitudinal/FirstLevel/SooEun/temp/NF102_1/Power_parameters.mat';
param = load(ParamTemplate);
roiMNI = param.parameters.rois.mni.coordinates;
NetsFile = '/net/parasite/HCP/Scripts/slab/PCA/FundDiffRepo/Data/Power_Nets.csv';
netscsv = dataset('File',NetsFile,'Delimiter',',');
netsmni = [netscsv.X netscsv.Y netscsv.Z netscsv.Network];
[a,b] = ismember(roiMNI,netsmni(:,1:3),'rows');
nets = netsmni(b,4);

Nets = max(nets);

for iNet = 1:(Nets-1)
    
    for jNet = (iNet+1):Netsak
        fprintf(1,'Nets %d,%d\n',iNet,jNet);
        
        netmask = zeros(Nodes,Nodes);
        netmask(nets==iNet,nets==iNet) = 1;
        netmask(nets==iNet,nets==jNet) = 1;
        netmask(nets==jNet,nets==iNet) = 1;
        netmask(nets==jNet,nets==jNet) = 1;
        netmask = mc_flatten_upper_triangle(netmask)>0;
        
        NumComp = 5;
        NumComps = size(NumComp,2);
        
        featuremat_abcd2 = abcd.featuremat_abcd(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_abcd2);
        Abig = score;
        icasig = coeff';

        %mean center 
        mu = mean(featuremat_abcd2);
        x = bsxfun(@minus,featuremat_abcd2,mu);

        %calculate expressions for each subject
        A = (pinv(icasig')*x')';
           
        comp = size(A,2);
        
        for iFold = 1:nFold
        fprintf(1,'.');
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
        fprintf(1,'\n');

        [x_pos_abcd,y_pos_abcd,t_pos_abcd,auc_pos_abcd(iNet,jNet)] = perfcurve(pheno, pheno_pred_abcd,2);
 
        featuremat_hcp2 = hcp.featuremat_hcp(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_hcp2);
        Abig = score;
        icasig = coeff';
        NumComp = 5;
        NumComps = size(NumComp,2);

        %mean center 
        mu = mean(featuremat_hcp2);
        x = bsxfun(@minus,featuremat_hcp2,mu);

        %calculate expressions for each subject
        A = (pinv(icasig')*x')';
        comp = size(A,2);

        pheno_pred_hcp = zeros(n,p,NumComps);

        for iFold = 1:nFold
        fprintf(1,'.');
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
        fprintf(1,'\n');

        [x_pos_hcp,y_pos_hcp,t_pos_hcp,auc_pos_hcp(iNet,jNet)] = perfcurve(pheno, pheno_pred_hcp,2);
        
        
        featuremat_2 = featuremat(:,netmask);
        [coeff, score, latent, ~, exp] = pca(featuremat_2);
        Abig = score;
        icasig = coeff';
       
        pheno_pred = zeros(n,p,NumComps);
        for iFold = 1:nFold
        fprintf(1,'.');
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
        fprintf(1,'\n');
        [x_pos,y_pos,t_pos,auc_pos(iNet,jNet)] = perfcurve(pheno, pheno_pred,2);
        figure; plot(x_pos,y_pos,'r-')

    end
    
end
