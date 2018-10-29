% Dissimilarity distance for all-vs-all viral proteins
% Part of: Dissimilarity-based Negative Sampling


%% STEP INPUT: Initialization
BLOUSUMno=30;
ResultFile=['Distance' int2str(BLOUSUMno) '.mat'];

%% STEP 1: Load Viral Protein Sequences
load Mentha_All_Vs.mat ProteinB_Unique
[~, SeqV]=fastaread('All_Viral_interacting.fasta');
VpCount=length(SeqV);

%% STEP 2: Global Alignment
ScoringMatrix=['BLOSUM' num2str(BLOUSUMno)];
ScoreVps=zeros(VpCount,VpCount);
parfor r=1:VpCount
    for c=1:VpCount
        ScoreVps(r,c)= nwalign(SeqV{1,r},SeqV{1,c},'ScoringMatrix',ScoringMatrix); %#ok<PFBNS>
    end
end

%% STEP 3: Remove Outliers
ScoreVps(:,133)=0;
ScoreVps(133,:)=0;
ScoreVps(:,333)=0;
ScoreVps(333,:)=0;

%% STEP4: Normalize for Each Viral Protein (row-wise)
Weights=zeros(VpCount,VpCount);
for r=1:VpCount
    Row=ScoreVps(r,:);
    %Normalize on [0-1] scale, then reverse it (1-x)
    Weights(r,:)=1-((Row-min(Row))/(max(Row)-min(Row)));
end

%% STEP OUTPUT: Save
save(ResultFile,'VpCount','BLOUSUMno','Weights')