% Cross Validation over the +/- interactions data set
% Runs over T, C, and Gamma


%% STEP INPUT: Initialization
clear
Threshold_T=0.9:-0.1:0 %#ok<NOPTS> % Dissimilarity threshold (T)
Kfolds=5 %#ok<NOPTS> %k-fold in CV
load ParIndex.mat
CGs=ParIndexCG2; % List of C and Gamma for the SVM model

%% Loop for each T value
for Fi=1:length(Threshold_T)
    
    %% STEP 1: Loop Inputs
    FeatureFile=['Features_C1_T' num2str(Threshold_T(Fi))  '.mat'] %#ok<NOPTS>
    ResultFile= ['Result_CV_T' num2str(Threshold_T(Fi))  '.mat'] %#ok<NOPTS>
    
    %% STEP 2: Load Features File
    load(FeatureFile,'TrainingSp','Training_LabelSp')
    SetCount=length(TrainingSp) %#ok<NOPTS>
    
    %% STEP 3: Prepare the Data Set
    Crnt_Train=[];
    Crnt_Label=[];
    for si=1:SetCount
        Crnt_Train=[Crnt_Train;TrainingSp{si,1}]; %#ok<AGROW>
        Crnt_Label=[Crnt_Label;Training_LabelSp{si,1}];  %#ok<AGROW>
    end
    
    %%  STEP 4: Standarization
    % Standarize Training and Testing
    
    FeatureCount=size(Crnt_Train,2);
    Mean_Train=zeros(1,FeatureCount);
    SD_Train=zeros(1,FeatureCount);
    % Mu and Segma of training
    Mean_Train(1,:)= mean(Crnt_Train,1);
    SD_Train(1,:)=std(Crnt_Train,0,1);
    CrntMean=Mean_Train(1,:);
    CrntSD=SD_Train(1,:);
    L_Tr=size(Crnt_Train,1);
    MuMatrix_Train=repmat(CrntMean,L_Tr,1);
    SDMatrix_Train=repmat(CrntSD,L_Tr,1);
    % Standarize
    Crnt_Train=(Crnt_Train-MuMatrix_Train)./SDMatrix_Train;
    %Remove NaNs (created by feature values of all 0s for a specific feature
    %acrros all feature vectors)
    Crnt_Train(find(isnan(Crnt_Train)==1))=0; %#ok<FNDSB>
    %Remove inf (created by scaling transfer to testing data)
    Crnt_Train(find(isinf(Crnt_Train)==1))=0; %#ok<FNDSB>
    
    
    %% STEP 5: Loop Initialization
    Acc_CV=zeros(size(CGs,1),1);
    CVTime=zeros(size(CGs,1),1);
    
    %% Looping Over the Model Paramters (here C, and Gamma)
    parfor BB=1:size(CGs,1)
        %% STEP 6: Set the Model Paramters
        CrntC=CGs(BB,1);
        CrntG=CGs(BB,2); %#ok<PFBNS>
        
        % LIBSVM SVM model paramter defintions:
        %{
-d degree : set degree in kernel function (default 3)
-g gamma : set gamma in kernel function (default 1/num_features)
-r coef0 : set coef0 in kernel function (default 0)
-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
-m cachesize : set cache memory size in MB (default 100)
-e epsilon : set tolerance of termination criterion (default 0.001)
-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)
-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)
-v n: n-fold cross validation mode
-q : quiet mode (no outputs)

        %}
        ParVect= [' -v '  int2str(Kfolds)  ... this is CV
            ' -c '  num2str(CrntC)...
            '  -t ' num2str(2) ...
            ' -s '  num2str(0) ...
            ' -h ' num2str(0) ...
            ' -m ' num2str(10000) ...
            ' -g ' num2str(CrntG) ...
            ]
        
        %% STEP 7: Do the Cross-Validation
        TotalTrainTime=tic;
        Acc_CV(BB,1)=svmtrainX(Crnt_Label,sparse(double(Crnt_Train)), ParVect); % Accuracy of CV
        CVTime(BB,1)=toc(TotalTrainTime) % Time of CV
    end
    
    %% STEP OUTPUT: Save
    save(ResultFile,'CVTime','Acc_CV','FeatureFile','Kfolds','CGs')
end
