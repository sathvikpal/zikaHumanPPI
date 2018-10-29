% Paritioning (Criterion 1) - Negative Sampling - Feature extraction

%% STEP INPUT: Initialization
global V_Seq
load FastaSeqData.mat V_Seq
Threshold_T=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
N2P_Ratio_train=1;% Negative:Positive Examples Ratio in Training
BLOSUMno=30;
%% Loop over T
for Fi=1:length(Threshold_T)
    %% STEP 1: Loop Paramters
    FileName=['Features_C1_T' num2str(CutList(Fi)) '.mat'];
    CutOFF=Threshold_T(Fi);
    %% STEP 2: Feature Extraction Paramters
    SamplingMode='CutRandom'; %Disimilarity-based Negative Sampling
    Description='Sampling randomly with Cutoff weight <T';
    Descriptor='Conjoint';
    if strcmp(Descriptor,'Conjoint')
        clear D_par
        D_par.Norm_Type='Shen'; % None, Shen, Cui
        D_par.AA= 'Shen';% Shen, Rand1 Rand2 Rand3 None
        D_par.kmer=3;%  % k [1:5], Larger k takes bulck of time and blocks memory
        D_par.Interface=0;
    end
    
    %% STEP 3:  Partitioning - Criterion 1
    % Group numbers in the comments are for organizing purpose only
    TrainSet=cell(0,0);
    TrainSet{end+1,1}=[134,166]; %#ok<*SAGROW> 
    TrainSet{end+1,1}=120; % GROUP 1: Paramyxo (1)
    TrainSet{end+1,1}=128; % GROUP 1: Paramyxo (62)
    TrainSet{end+1,1}=67; % GROUP 1: Paramyxo (12)
    TrainSet{end+1,1}=66; % GROUP 1: Paramyxo (4)
    TrainSet{end+1,1}=68; % GROUP 2: Penuma(222) 
    TrainSet{end+1,1}=[130,135]; % GROUP 3: Filo(90,1) Ebola
    TrainSet{end+1,1}=[111]; % GROUP 3: Filo(23)Marburg
    TrainSet{end+1,1}=[169]; % GROUP 4a: OrthoBunya(83)
    TrainSet{end+1,1}=[73]; % GROUP 4b: Thlebo(26)
    TrainSet{end+1,1}=[101]; % GROUP 4b: Thlebo(50)
    TrainSet{end+1,1}=[106,54]; % GROUP 5: Flavi(107,1)
    TrainSet{end+1,1}=[55]; % GROUP 5: Flavi(1) 
    TrainSet{end+1,1}=[62,148]; % GROUP 6:HCV-subtype (4,93)
    TrainSet{end+1,1}=[59,61]; % GROUP 6:HCV-subtype (2,38)
    TrainSet{end+1,1}=58; % GROUP 6:HCV-subtype (3)
    TrainSet{end+1,1}=[60,63,107,143,158,159]; % GROUP 6:HCV-subtype(1,2,1,28,4,6)
    TrainSet{end+1,1}=[26,100]; % GROUP 7: Adeno(18,46)_
    TrainSet{end+1,1}=[27,116]; % GROUP 7: Adeno (1,4)
    TrainSet{end+1,1}=[131]; % GROUP 7: Adeno (1)
    TrainSet{end+1,1}=[28 ]; % GROUP 7: Adeno(1)
    TrainSet{end+1,1}=[29,99]; % GROUP 7: Adeno (1,16)
    TrainSet{end+1,1}=[150,152,161,153,154,165,126,136];% GROUP 8: flu A H1N1(61,6,1,3,70,3,23,284)
    TrainSet{end+1,1}=[132]; % GROUP 8: flu A H5N1 Hongkong (108)
    TrainSet{end+1,1}= [151,157,155,162]; % GROUP 8: flu A H3N2 (101,1,2,1)
    TrainSet{end+1,1}=[104,2,3,1]; % GROUP 9: Chordopox(7,9,156,6)
    TrainSet{end+1,1}=[4];  % GROUP 9: Chordopoxng (16)
    TrainSet{end+1,1}=[41];% GROUP 10: Pabilloma-Alpha (2)
    TrainSet{end+1,1}=[38,145]; % GROUP 10: Pabilloma-Alpha (1,45)
    TrainSet{end+1,1}=[39,37,144];% GROUP 10: Pabilloma-Alpha  (1,12,114)
    TrainSet{end+1,1}=[35,40];% GROUP 10: Pabilloma-Alpha (36,16)
    TrainSet{end+1,1}=[42];% GROUP 10: Pabilloma-Kabba  (3)
    TrainSet{end+1,1}=[34];% GROUP 10: Pabilloma-Beta  (8)
    TrainSet{end+1,1}=[36];% GROUP 10: Pabilloma-Mu  (7)
    %Group 11: Herpses Beta
    TrainSet{end+1,1}=[18];% (2)
    TrainSet{end+1,1}=[19];%(1)
    TrainSet{end+1,1}=[14,15,16];%(62,2,19)
    % GROUP 12: Herpes Gamma
    TrainSet{end+1,1}=[21,124,20];%(103.13,297)
    TrainSet{end+1,1}=[115,160,170];%(20,115,10)
    % GROUP 13: Herpes Alpha 
    TrainSet{end+1,1}=[5,6,7,8,9];% (120,225,1,1,)
    TrainSet{end+1,1}=[10,11];%(2,7)
    % GROUP 14': HIV 1 
    TrainSet{end+1,1}=[149,86,84,82,80,77]; %(3,2,80,1,4,39)
    TrainSet{end+1,1}=[79,83,85,142] ;% (5,4,96,4)
    TrainSet{end+1,1}=[76] ;%(127)
    TrainSet{end+1,1}=[156];%#ok<*NBRAK> % (5)
    TrainSet{end+1,1}=[81];% (3)
    TrainSet{end+1,1}=[78];% (14)
    TrainSet{end+1,1}=[87]; % Group 15:HIV2 (8)
    TrainSet{end+1,1}=[88]; % Group 15:HIV2 (4)
    
    %% STEP 4: Converstion from Subsets TaxonIndex in TrainSet to TaxonID
    load Mentha_All_Vs.mat VirusID_Maper
    for t1=1:length(TrainSet)
        Crnt=TrainSet{t1,1};
        for c=1:length(Crnt)
            Crnt(c)=VirusID_Maper(Crnt(c));
        end
        TrainSet{t1,1}=Crnt;
    end
    clear VirusID_Maper Crnt    
    
    %% STEP 5: Negative Sampling
    % Data
    global Interactions_Table H_Hdr H_Seq V_Hdr %#ok<*TLEV>
    load Mentha_All_Vs.mat   Sub_Unique Tax_Unique  Interactions_Table
    load FastaSeqData.mat H_Hdr H_Seq  V_Hdr
    % Samling Modes
    %% Mode 1: Separation by Subset
    if ~strcmp(SamplingMode,'Proportional') && ~strcmp(SamplingMode,'CutRandom')
        NowIs='generating negative sample'; %#ok<*NOPTS>
        % Separate before sampling
        [P_perV_List_Train,H_NonInter_Train]= Negative_Sampling_BySet(SamplingMode,TrainSet); % V_List maps to ProteinB_Unique    
        [H_Inter_Seq_Train,H_NonInter_Seq_Train,~]=PrepareSetsSeq_BySet(P_perV_List_Train,H_NonInter_Train,N2P_Ratio_train);
    
    %% Mode 2: % Proportional Random Sampling
    elseif  strcmp(SamplingMode,'Proportional')           
        [P_perV_List_Train,H_NonInter_Train,H_Inter_Seq_Train,H_NonInter_Seq_Train]= ...
            Negative_Sampling_Proportional(TrainSet,N2P_Ratio_train,BLOSUMno,CutOFF); % V_List maps to ProteinB_Unique
        
    %% Mode 3: Disimilarity-based Negative Sampling
    elseif strcmp(SamplingMode,'CutRandom') 
        [P_perV_List_Train,H_NonInter_Train,H_Inter_Seq_Train,H_NonInter_Seq_Train]= ...
            Negative_Sampling_mod(TrainSet,N2P_Ratio_train,BLOSUMno,CutOFF); % V_List maps to ProteinB_Unique  
    end
   
    %% STEP 6: Feature Extraction  
    % Initialization
    FeatureVectors_Train_P=[];
    FeatureVectors_Train_N=[];
    IDs_Train_P=[];
    IDs_Train_N=[];
    SetCount=size(P_perV_List_Train,1);
    Vcount=size(P_perV_List_Train,2);    
    Training_LabelSp=cell(0,0);
    TrainingSp=cell(0,0);
    
    %% Loop over the Subsets
    parfor s=1:SetCount 
        
        %% Generating +ve Interaction Feature Vectors and Labels vector
        F_Vect_Temp=[];
        ID_Vect_Temp=[];
        for v=1:Vcount
            CrntV_Vps=P_perV_List_Train{s,v};
            for p=1:length(CrntV_Vps)
                Crnt_Vp=CrntV_Vps(p);
                Vp_H_Seq=H_Inter_Seq_Train{s,v,p};
                Crnt_F= single(FeatureExtraction_Mentha_SVM_ByVp(Crnt_Vp,Vp_H_Seq,Descriptor,D_par));
                F_Vect_Temp=[F_Vect_Temp;Crnt_F];
                Crnt_IDs=[repmat(Crnt_Vp,[size(Vp_H_Seq,1),1])];
                ID_Vect_Temp=[ID_Vect_Temp;Crnt_IDs];
            end
        end
        FeatureVectors_Train_P=F_Vect_Temp;
        IDs_Train_P=ID_Vect_Temp;
        
        %% Generating -ve Interaction Feature Vectors and Labels vector
        TrainNegTime=tic;
        F_Vect_Temp=[];
        for v=1:Vcount
            CrntV_Vps=P_perV_List_Train{s,v};
            for p=1:length(CrntV_Vps)
                Crnt_Vp=CrntV_Vps(p);
                Vp_H_Seq=H_NonInter_Seq_Train{s,v,p};
                Crnt_F= single(FeatureExtraction_Mentha_SVM_ByVp(Crnt_Vp,Vp_H_Seq,Descriptor,D_par));
                F_Vect_Temp=[F_Vect_Temp;Crnt_F];
            end
        end
        FeatureVectors_Train_N=F_Vect_Temp;
        TrainNegTime=toc(TrainNegTime);
      
        %% Combine + and - for each subset     
        Training_LabelSp{s,1}=[ones(size(FeatureVectors_Train_P,1),1);-1*ones(size(FeatureVectors_Train_N,1),1)];
        TrainingSp{s,1}=[FeatureVectors_Train_P;FeatureVectors_Train_N]; 
   
    end % End of loop for a subset
    
    %% STEP OUTPUT: Save
    save(FileName,'CutOFF','Description', 'TrainingSp','Training_LabelSp')     
end