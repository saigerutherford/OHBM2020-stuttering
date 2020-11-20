%% Add needed paths for loading connectomes & phenotype data

addpath('/usr/local/MethodsCore/matlabScripts/')
addpath('/usr/local/MethodsCore/matlabScripts/Takgraph/')

ParamTemplate = '/net/pepper/youth_longitudinal/FirstLevel/SooEun/temp/NF102_1/Power_parameters.mat';
CorrTemplate = '/net/pepper/youth_longitudinal/FirstLevel/SooEun/temp/[Subject]/Power_corr.mat';

% Path to phenotype file
PhenotypeFile = '/net/pepper/youth_longitudinal/Scripts/SooEun/LogReg/sooeunscanfile_visit1included_allcov.csv'; % NvsS1

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

% %% Set variables for consensus map visualization
NetsFile = '/net/parasite/HCP/Scripts/slab/PCA/FundDiffRepo/Data/Power_Nets.csv';

netscsv = dataset('File',NetsFile,'Delimiter',',');
netsmni = [netscsv.X netscsv.Y netscsv.Z netscsv.Network];
[a,b] = ismember(roiMNI,netsmni(:,1:3),'rows');
nets = netsmni(b,4);

%% Clear unused variables and save workspace
clear a b bad CorrTemplate dat Edges iN iSubject n_train n_test nan netscsv NetsFile netsmni Nodes p param ParamPath ...
    ParamPathCheck ParamTemplate path PhenotypeFile r replace roiMNI subs z

save('InSample.mat')