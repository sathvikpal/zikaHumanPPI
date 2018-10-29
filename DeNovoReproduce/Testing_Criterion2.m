% Training - Testing on Criterion 1 Partitioning

%% STEP INPUT: Testing Paramters
clear
CGs=[10,0.001]; % Optimal paramters found in the CV
Threshold_T=[0:0.1:0.9]; %#ok<*NBRAK> % Dissimilarity threshold

for Cu=1:length(Threshold_T)
     %% STEP 1: Loop Paramters
    CrntCut=Threshold_T(Cu);
    FeatureFile=['Features_C2_T' num2str(CrntCut) '.mat'];
    ResultFile =['Result_C2_Prob_MOD_T' num2str(CrntCut) '.mat'];
    
    %% STEP 2: Load Loop Data
    load(FeatureFile,'Description','TrainingSp','Training_LabelSp')
        
    %% STEP 3: Initialization
    SetCount=length(TrainingSp);
    ParVect_cell=cell(size(CGs,1),1);
    endtime_cell=zeros(size(CGs,1),1);
    Stat_Ts_cell=cell(size(CGs,1),1);
    Averages_cell=zeros(size(CGs,1),7);
    Stat_Tr_cell=cell(size(CGs,1),1);
    
    %% Loop over the input paramter pairs (C,Gamma)
    for CrntPointer=1:size(CGs,1)
        StartTime=tic;
        CrntC=CGs(CrntPointer,1);
        CrntG=CGs(CrntPointer,2);
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
        ParVect= [' -c '  num2str(CrntC)...
            '  -t ' num2str(2) ...
            ' -s '  num2str(0) ...
            ' -h ' num2str(0) ...
            ' -m ' num2str(10000) ...
            ' -g ' num2str(CrntG) ...
            ' -b ' num2str(1) ...
            ];
        
       
        %% STEP 4: Training
        %% STEP 4.1 Initializations
        Tr_Sens=zeros(SetCount,1);
        Tr_Spec=zeros(SetCount,1);
        Tr_Acc=zeros(SetCount,1);
        Tr_MCC=zeros(SetCount,1);
        Tr_LibAcc=zeros(SetCount,1);
        Tr_Size=zeros(SetCount,1);
        Ts_Sens=zeros(SetCount,1);
        Ts_Spec=zeros(SetCount,1);
        Ts_Acc=zeros(SetCount,1);
        Ts_MCC=zeros(SetCount,1);
        Ts_LibAcc=zeros(SetCount,1);
        Ts_Size=zeros(SetCount,1);
        SVratio=zeros(SetCount,1);
        ConfideneceMat=cell(SetCount,1);
        %% Loop over each data subset
        parfor s=1:SetCount
            %% STEP 4.2: Define Training and Testing
            % leave-one-out 
             Included=setdiff([1:SetCount],s);
            Crnt_Train=[];
            Crnt_Label=[];
            for si=1:SetCount-1
                Crnt_Train=[Crnt_Train;TrainingSp{Included(si),1}]; %#ok<*PFBNS>
                Crnt_Label=[Crnt_Label;Training_LabelSp{Included(si),1}];
            end
            CrntTest=TrainingSp{s,1};
            Testing_Label=Training_LabelSp{s,1};
            
            
            %%  STEP 4.3: Standarization
            FeatureCount=size(Crnt_Train,2);
            Mean_Train=zeros(1,FeatureCount);
            SD_Train=zeros(1,FeatureCount);            
            % Mu and Segma            
            Mean_Train(1,:)= mean(Crnt_Train,1);
            SD_Train(1,:)=std(Crnt_Train,0,1);            
            CrntMean=Mean_Train(1,:);
            CrntSD=SD_Train(1,:);
            L_Tr=size(Crnt_Train,1);
            MuMatrix_Train=repmat(CrntMean,L_Tr,1);
            SDMatrix_Train=repmat(CrntSD,L_Tr,1);
            L_Ts=size(CrntTest,1);
            MuMatrix_Test=repmat(CrntMean,L_Ts,1);
            SDMatrix_Test=repmat(CrntSD,L_Ts,1);
            % Standarize BOTH Training and Testing
            Crnt_Train=(Crnt_Train-MuMatrix_Train)./SDMatrix_Train;
            CrntTest=(CrntTest-MuMatrix_Test)./SDMatrix_Test;
            %Remove NaNs
            Crnt_Train(find(isnan(Crnt_Train)==1))=0; %#ok<*FNDSB>
            CrntTest(find(isnan(CrntTest)==1))=0;
            %Remove inf % created by scaling transfer to testing data
            Crnt_Train(find(isinf(Crnt_Train)==1))=0;
            CrntTest(find(isinf(CrntTest)==1))=0;         
            
            %% STEP 4.4 Training SVM
            TotalTrainTime=tic;
            ModelTemp=svmtrainX(Crnt_Label,sparse(double(Crnt_Train)), ParVect);
            SVratio(s,1)=ModelTemp.totalSV/length(Crnt_Label);
            TotalTrainTime=toc(TotalTrainTime);
            
            %% STEP 4.5: Training Accuracy            
            [predicted_label,Tr_Accuracy,~] = svmpredict(Crnt_Label, sparse(double(Crnt_Train)), ModelTemp,'-b 1');
            CM=confusionmat(Crnt_Label,predicted_label,'order',[1 -1]);
            TP=CM(1,1);
            TN=CM(2,2);
            FN=CM(1,2);
            FP=CM(2,1);
            Tr_Sens(s,1)=TP/(TP+FN);
            Tr_Spec(s,1)=TN/(FP+TN);
            Tr_Acc(s,1)=(TP+TN)/sum(sum(CM));
            Tr_MCC(s,1)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
            Tr_LibAcc(s,1)=Tr_Accuracy(1);
            Tr_Size(s,1)=length(Crnt_Label);
             
            %% STEP 5: Testing and Testing  Accuracy
            [predicted_label,Ts_Accuracy,ConfideneceMat{s}] = svmpredict(Testing_Label, sparse(double(CrntTest)), ModelTemp,'-b 1');
            CM=confusionmat(Testing_Label,predicted_label,'order',[1 -1]);
            TP=CM(1,1);
            TN=CM(2,2);
            FN=CM(1,2);
            FP=CM(2,1);
            Ts_Sens(s,1)=TP/(TP+FN);
            Ts_Spec(s,1)=TN/(FP+TN);
            Ts_Acc(s,1)=(TP+TN)/sum(sum(CM));
            Ts_MCC(s,1)=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
            Ts_Size(s,1)=length(Testing_Label);
            
        end
        
        %% STEP 6: Prepare Results
        Stat_Tr=[Tr_Acc,Tr_Sens,Tr_Spec,Tr_MCC,Tr_Size,SVratio];
        Stat_Ts=[Ts_Acc,Ts_Sens,Ts_Spec,Ts_MCC,Ts_Size]; 
        Averages(1)=sum(Ts_Acc.*Ts_Size)/sum(Ts_Size); % Weighted accuracy
        Averages(2)= sum(Ts_Sens.*Ts_Size)/sum(Ts_Size); % Weighted sensitivity
        Averages(3)= sum(Ts_Sens)/length(Ts_Size); % Avergae sensitivity
        Averages(4)= sum(Ts_Spec.*Ts_Size)/sum(Ts_Size); % Weighted sepecificity
        Averages(5)= sum(Ts_Spec)/length(Ts_Size); % Average sepecificity
        Averages(6)=mean(Ts_Acc); % Average accuracy
        Averages(7)=sum(SVratio.*Ts_Size)/sum(Ts_Size);
        ParVect_cell{CrntPointer,1}=ParVect;
        endtime_cell(CrntPointer,1)=toc(StartTime);
        Stat_Ts_cell{CrntPointer,1}=Stat_Ts;
        Averages_cell(CrntPointer,:)=Averages;
        Stat_Tr_cell{CrntPointer,1}=Stat_Tr;
        
    end
    
    %% STEP OUTPUT: Save
    save(ResultFile,'ConfideneceMat','FeatureFile','ParVect_cell','Averages_cell','Stat_Ts_cell','Stat_Tr_cell','endtime_cell') %ConcatResult % Best Par =10,10: 72,61,84
end
