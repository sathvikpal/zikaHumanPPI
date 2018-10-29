% Paritioning (Criterion 1) - Negative Sampling - Feature extraction

%% STEP INPUT: Initialization
global V_Seq
load FastaSeqData.mat V_Seq
Threshold_T=[0:0.1:0.9]; %#ok<*NBRAK>
N2P_Ratio_train=1;% Negative:Positive Examples Ratio in Training
BLOSUMno=30;
%% Loop over T
for Fi=1:length(Threshold_T)
     %% STEP 1: Loop Paramters
    FileName=['Features_C2_T_' num2str(CutList(Fi)) '.mat']
    CutOFF=Threshold_T(Fi);
    Description=['A family subset, Sampling randomly with Cutoff weight <' num2str(Threshold_T(Fi))];
     
    %% STEP 2: Feature Extraction Paramters
    Descriptor='Conjoint';
    SamplingMode='CutRandom';   
    if strcmp(Descriptor,'Conjoint')
        clear D_par
        D_par.Norm_Type='Shen'; % None, Shen, Cui
        D_par.AA= 'Shen';% Shen, Rand1 Rand2 Rand3 None
        D_par.kmer=3;%   % k [1:5], Larger K takes time on laptop abd blocks memory
        D_par.Interface=0;
    end
    
    
    %% STEP 3:  Partitioning - Criterion 1
    % Group nunmbers in the comments are for organizing purpose only
    TrainSet=cell(0,0);
    % Adding species separately
    % Family 1: Paramyxo (448,13,1,62,12,4,222)
    TrainSet{end+1,1}=[134,166,120,128,67,66,68]; %#ok<*SAGROW> 
    % Family 2: Filo(90,1,23)Ebola & Murburg
    TrainSet{end+1,1}=[130,135,111];
    % Family 3: Bunya (83,26,50)
    TrainSet{end+1,1}=[169,73,101];   
    % Family 4: Flavi (107,1,1,4,93,2,38,3,1,2,1,28,4,6
    TrainSet{end+1,1}=[106,54,55,62,148,59,61,58,60,63,107,143,158,159];   
    % Family 5: Adeno (18,46,1,4,1,1,1,16)
    TrainSet{end+1,1}=[26,100,27,116,131,28,29,99 ]; % GROUP 7: Adeno(18,46)_ 
    % Family 6: Orthomyxo H1N1(61,6,1,3,70,3,23,284)H5N1(108)H3N2 (101,1,2,1)
    TrainSet{end+1,1}=[150,152,161,153,154,165,126,136,132,151,157,155,162];% GROUP 8: flu A  
    % Family 7: Chordpox (7,9,156,6,16)
    TrainSet{end+1,1}=[104,2,3,1,4];
    % Family 8: Papiloma  (2,1,45,1,12,114,46,16,3,8,7)
    TrainSet{end+1,1}=[41,38,145,39,37,144,35,40,42,34,36];% GROUP 10: Pabilloma-Alpha (2)
    % Family 9: Herpses (2,1,62,2,19,103,13,297,20,115,10,120,225,1,1,2,7)
    TrainSet{end+1,1}=[18,19,14,15,16,21,124,20,115,160,170,5,6,7,8,9,10,11];% (2)
    % Family 10: Retro (3,2,80,1,4,39,5,4,96,4,127,5,3,14,8,4)
    TrainSet{end+1,1}=[149,86,84,82,80,77,79,83,85,142,76,156,81,78,87,88];
    
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
    global Interactions_Table H_Hdr H_Seq V_Hdr %#ok<TLEV>
    load Mentha_All_Vs.mat   Sub_Unique Tax_Unique  Interactions_Table
    load FastaSeqData.mat H_Hdr H_Seq  V_Hdr
    % Samling Modes
    %% Mode 1: Separation by Subset
    if ~strcmp(SamplingMode,'Proportional') && ~strcmp(SamplingMode,'CutRandom')
        % Separate Before Sampling
        [P_perV_List_Train,H_NonInter_Train]= Negative_Sampling_BySet(SamplingMode,TrainSet); % V_List maps to ProteinB_Unique
        % Prepare Sets of Sequences
        [H_Inter_Seq_Train,H_NonInter_Seq_Train,~]=PrepareSetsSeq_BySet(P_perV_List_Train,H_NonInter_Train,N2P_Ratio_train);
    
    %% Mode 2: Proportional Random Sampling     
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
    SetCount= size(P_perV_List_Train,1);
    Vcount=size(P_perV_List_Train,2);
    Training_LabelSp=cell(0,0);
    TrainingSp=cell(0,0);
    
    %% Loop over the Subsets
    parfor s=1:SetCount % parfor
        
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
        TrainNegTime=toc(TrainNegTime)

           %% Combine + and - for each subset         
        Training_LabelSp{s,1}=[ones(size(FeatureVectors_Train_P,1),1);-1*ones(size(FeatureVectors_Train_N,1),1)];
        TrainingSp{s,1}=[FeatureVectors_Train_P;FeatureVectors_Train_N];
        
    end % End of loop for a subset
    
    %% STEP OUTPUT: Save
    save(FileName,'CutOFF','Description', 'TrainingSp','Training_LabelSp')
end