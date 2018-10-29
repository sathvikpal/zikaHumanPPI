function [P_perV_List,H_NonInter,H_Inter_Seq,H_NonInter_Seq] = ...
    Negative_Sampling_Proportional(Vset,N2P_Ratio,BLOUSUMno,CutOFF)
% Dissimilarity-based Negative Sampling with the outliers masked in the input
%{
Input:
    Vset: cell(n,1) for n different Array sets, to process each individually
     Vset elements MUST be indexed to Tax_Unique
  
Output:
    P_perV_List: cell(n,max(#V in Sets))-> Array (Vp IDs for each V set)
                    row: Vset#, col:V#
    H_NonInter: H_NonInter{Snum,Vnum,p}=NegSet_P
%}


% Generate sets of non-interacting proteins wwith each virus/viral protein
%{
Input:
    SamplingMode: How to generate negative set
    ['Random','RandomNonSimilar']
    Vset: cell(n,1) for n different Array sets, to process each individually
     Vset elements MUST be indexed to Tax_Unique
  
Output:
    P_perV_List: cell(n,max(#V in Sets))-> Array (Vp IDs for each V set)
                    row: Vset#, col:V#
    H_NonInter: H_NonInter{Snum,Vnum,p}=NegSet_P

%}
%% Data

% Interactions
global Interactions_Table
% Interaction Table: rows (interactions) columns (H|V)
Hcol=Interactions_Table(:,1);  % Human proteins indexed to ProteinA_Unique
Vpcol=Interactions_Table(:,2); % Viral Protiens IDs indexed to ProteinB_Unique
Vcol=Interactions_Table(:,3);  % Virus IDs from Tax_Unique = TaxonB

% Dissimilarity Distances
WeightsFile=['Distance' int2str(BLOUSUMno) '.mat'];
load(WeightsFile,'VpCount','BLOUSUMno','Weights')
% Unique Vs arranged in order consistent with Protein(A,B)_Unique Indexing Indexing

% Initializations
P_perV_List=cell(0,0);
H_NonInter=cell(0,0,0);

%% Sampling
for Snum=1:length(Vset)
    % For each virus
    CrntSet=Vset{Snum,1};
    for Vnum=1:length(CrntSet)
        % List of Proteins in Each Viruses
        CrntV_Ps_Inx=find(Vcol==CrntSet(Vnum));
        CrntV_Vps_Pool=unique(Vpcol(CrntV_Ps_Inx)); %#ok<*FNDSB>
        P_perV_List{Snum,Vnum}=CrntV_Vps_Pool;
        % For each protein,
        for p=1:length(CrntV_Vps_Pool)
            CrntP=CrntV_Vps_Pool(p);
            %  Find its interactors
            HsofP_Inx=find(Vpcol==CrntP  & Vcol==CrntSet(Vnum));
            HsofP=Hcol(HsofP_Inx);
            % Generate -ve Set: Select randomly from ~Hs(P) after T
            Pos_count=length(HsofP);
            Neg_count=N2P_Ratio*Pos_count;
            NegSet_P=[];
            while length(setdiff(NegSet_P,HsofP))<Neg_count
                % Pick Vp' randomly with weights
                PropWeight=Weights(CrntP,:); %#ok<*NODEF>
                PropWeight(PropWeight<CutOFF)=0;
                Non_Vp= datasample([1:VpCount],1,'Weights',PropWeight); %#ok<*NBRAK>
                % Pick H from Vp' interactors uniformally
                Non_Vp_Hs=Hcol(find(Vpcol==Non_Vp));
                AddedH=Non_Vp_Hs(randi(length(Non_Vp_Hs),1));
                % Add H to non-Interacting for Vp
                NegSet_P=[NegSet_P;AddedH]; %#ok<*AGROW>
            end
            NegSet_P=setdiff(NegSet_P,HsofP); % Gives unique value
            % Add to the return H_NonInter
            H_NonInter{Snum,Vnum,p}=NegSet_P;
        end
    end
end

%% Prepare Sequences
P_List=P_perV_List;
% Read Fasta
global H_Hdr H_Seq  V_Hdr
% Indexing Fasta
% Uniport is 6 char in Header(4:9)
H_Index=cell(length(H_Hdr),1);
for h=1:length(H_Hdr)
    H_Index{h,1}=H_Hdr{1,h}(4:9);
end
V_Index=cell(length(V_Hdr),1);
for v=1:length(V_Hdr)
    V_Index{v,1}=V_Hdr{1,v}(4:9);
end
% Fatsa are sorted exactly as ProteinA_Unique and ProteinB_Unique

%% Positive Sequences
load Mentha_All_Vs.mat Interactions_Table % H|V
Hcol=Interactions_Table(:,1);
Vpcol=Interactions_Table(:,2);
SetsCount=size(P_List,1);
Vcount=size(P_List,2);
H_Inter_Seq=cell(0,0);
for S=1:SetsCount
    for V=1:Vcount
        % Viral proteins of Crnt Virus
        CrntV_ps=P_List{S,V};
        VpCount=length(CrntV_ps);
        for P=1:VpCount
            % For each viral protein in this virus
            CrntVp=CrntV_ps(P);
            % find its Human interactors
            CrntVp_Pos=Hcol(find(Vpcol==CrntVp));
            if ~isempty(CrntVp_Pos)
                TempCell=cell(length(CrntVp_Pos),1); % Cells inside Cells
                for hin=1:length(CrntVp_Pos)
                    % Add the sequnce of each human interactor for this Vp
                    TempCell{hin}=H_Seq{1,CrntVp_Pos(hin)};
                end
                H_Inter_Seq{S,V,P}=TempCell;
            end
        end
    end
end

%% Negative Sequences
% H NONinteractng Seq [Stored as Cells inside Cells ]
H_NonInter_Seq=cell(0,0);%cell(SetsCount,Vcount,VpCount);
H_NonInter_Mod=H_NonInter; %#ok<*NASGU>
for S=1:SetsCount
    for V=1:Vcount
        CrntV_ps=P_List{S,V};
        VpCount=length(CrntV_ps);
        for P=1:VpCount
            % The negative H set for this Vp
            CrntVp_Neg=H_NonInter{S,V,P};
            % Now prepare the reduced -ve set sequences
            if ~isempty(CrntVp_Neg)
                TempCell=cell(length(CrntVp_Neg),1); % Cells inside Cells
                for hin=1:length(CrntVp_Neg)
                    TempCell{hin,1}=H_Seq{1,CrntVp_Neg(hin)};
                end
                H_NonInter_Seq{S,V,P}=TempCell;
            end
        end
    end
end
end

