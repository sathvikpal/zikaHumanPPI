function FeatureVectors= FeatureExtraction_Mentha_SVM_ByVp(Vp,Vp_H_Seq,Descriptor,D_par)
% Dissimilarity-based Negative Sampling for a single viral protein
%{
INPUT:
    Vp: the viral protein to couple with its list of +/-interactors
    Vp_H_Seq{#,1}: Cell of Character Array of sequences for Hs interacting with Vp
    Descriptor,D_par: Defines the feature extraction scheme
OUTPUT:
    FeatureVectors(#,AAg^kmer) Matrix (kmer: number of AA to cinsider in a
    single window when calculating frequencies)
%}

% Data and Initializations
global V_Seq
% AA Definitions
CiuAA={'IVLM', 'FYW', 'HKR', 'DE', 'QNTP','ACGS'};
ShenAA={'AVG', 'ILFP', 'YMTS', 'HNQW', 'RK','DE','C'};
Rand1AA={'ACDE','FGHI','KLMN','PQRS','TVWY'};
Rand2AA={'ACQRM','FTVHG','KLID','PSEWNY'};
Rand3AA={'AC','FRMG','KLD','PSEY','I','QWN','TVH'};
NoAA={'A','C','F','R','M','G','K','L','D','P','S','E','Y','I','Q','W','N','T','V','H'};

% Account for functional similarities between two discriptors
if strcmp(Descriptor,'Interface')
    Descriptor='Conjoint';
end

% The study Main Discriptor:
if strcmp(Descriptor,'Conjoint')
    Norm_Type=D_par.Norm_Type;
    AA=D_par.AA;
    kmer=D_par.kmer;
    
    % Set AA-related paramters
    if strcmp(AA,'Cui')
        AA_class=CiuAA;
        AA_class_n=6;
    elseif strcmp(AA,'Shen')
        AA_class=ShenAA;
        AA_class_n=7;
    elseif strcmp(AA,'Rand1')
        AA_class=Rand1AA;
        AA_class_n=5;
    elseif strcmp(AA,'Rand2')
        AA_class=Rand2AA;
        AA_class_n=4;
    elseif strcmp(AA,'Rand3')
        AA_class=Rand3AA;
        AA_class_n=7;
    elseif strcmp(AA,'None')
        AA_class=NoAA;
        AA_class_n=20;
    end
    
    %% AA Clustering [map AA into one of the 7/6/... AA classes
    
    % Map the single Virus Protein given
    Vp_mapedSeq=[]; %#ok<*NASGU>
    CrntVpSeq=V_Seq{1,Vp};
    CrntVp_mapedSeq=zeros(length(CrntVpSeq),1);
    for n=1:length(CrntVpSeq)
        for AAc=1:AA_class_n
            if ~isempty(strfind(AA_class{AAc},CrntVpSeq(n)))
                CrntVp_mapedSeq(n)=AAc;
            end
        end
    end
    Vp_mapedSeq=CrntVp_mapedSeq;
    
    % Map Human Proteins Given
    h_count=size(Vp_H_Seq,1);
    H_mapedSeq=cell(h_count,1);
    for hin=1:h_count
        CrntSeq=Vp_H_Seq{hin,1}; % Character Array
        TempCell=zeros(1,length(CrntSeq));
        for n=1:length(CrntSeq)
            for AAc=1:AA_class_n
                SSS=CrntSeq(n);
                try
                    if ~isempty(strfind(AA_class{AAc},SSS))
                        TempCell(1,n)=AAc;
                    end
                catch
                    SSS
                end
            end
        end
        H_mapedSeq{hin,1}=TempCell;
    end
    
    %% Building Feature Vectors:
    %   Mapping to fixed-length feature vector (AA_classes^kmer), and Counting frequencies of k-mer windows
    
    % Generating All Kmers combinations
    Kmers_comb=cell(AA_class_n^kmer,1);
    kmer_Mat=unique(combnk(repmat([1:AA_class_n]',kmer,1),kmer),'rows');
    if kmer==1
        kmer_Mat=kmer_Mat';
    end
    for Trip=1:size(Kmers_comb,1)
        for k=1:kmer
            Kmers_comb{Trip,1}(k,1)=int2str(kmer_Mat(Trip,k));
        end
        Kmers_comb{Trip,1}=Kmers_comb{Trip,1}';
    end
    
    % Feature Vectors: Single Virus Protein given
    F_Vseq=zeros(1,AA_class_n^kmer);
    CrntVp_mapedSeq=Vp_mapedSeq;
    String_Vseq=(int2str(CrntVp_mapedSeq))';
    for t=1:length(Kmers_comb)% =AA_class_n^kmer
        F_Vseq(1,t)=length(strfind(String_Vseq,Kmers_comb{t,1}));
    end
    
    % Feature Vectors: Associated Human Proteins
    F_HSeq=zeros(h_count,AA_class_n^kmer);
    for hin=1:h_count
        CrntVH=H_mapedSeq{hin,1};
        TempStr=int2str(CrntVH');
        for K=1:length(Kmers_comb)
            F_HSeq(hin,K)=length(strfind(TempStr',Kmers_comb{K,1}));
        end
    end
    
    %% Normalization
    % Ciu2012: di={fi-min(f)/e^(max(f)-min(f))}  -1
    % Shen2007: di=fi-min(f)/(max(f)-min(f))
    
    % Virus Protein
    F_Vseq_norm=NormalizeSVM(Norm_Type,F_Vseq);
    % H_Proteins
    F_HSeq_norm=zeros(size(F_HSeq,1),size(F_HSeq,2));
    for h=1:size(F_HSeq,1)
        F_HSeq_norm(h,:)=NormalizeSVM(Norm_Type,F_HSeq(h,:));
    end
    
    % Combine feature Vectors (Vp,Hp)
    % SVM data entry: row = sample, col= feature
    
    %% Descriptor 1: Interface
    if D_par.Interface==0 %normal Dijoint, no interface
        FeatureVectors=[repmat(F_Vseq_norm,size(F_HSeq,1),1),F_HSeq_norm];
    elseif D_par.Interface==1 %Interface Discriptor
        FeatureVectors=zeros(size(F_HSeq,1),length(F_Vseq_norm)^2);
        Utrans=F_Vseq_norm';
        for rrr=1:size(F_HSeq,1) %Go through every pair
            UV_Mat=(Utrans*F_HSeq_norm(rrr,:))';
            %convert to single to reduce memory usage
            FeatureVectors(rrr,:)=single(UV_Mat(:));
        end
    end
end

%% Descriptor 2: AAcomposition
if strcmp(Descriptor,'AAcomposition')
    
    % Frequency of each AA in each protein  AND Normalize
    AAs='GPAVLIMCFYWHKRQNEDST';
    % The single Virus Protein given
    Vp_mapedSeq=[];
    CrntVpSeq=V_Seq{1,Vp}; %Vp is a functin input of virus protein Index
    CrntVp_mapedSeq=zeros(1,length(AAs));
    for n=1:length(AAs)
        CrntVp_mapedSeq(1,n)=(length(strfind(CrntVpSeq,AAs(n))))/length(CrntVpSeq);
    end
    % Assigning variable names to match previously built code
    Vp_mapedSeq=CrntVp_mapedSeq;
    F_Vseq=Vp_mapedSeq;
    F_Vseq_norm=F_Vseq;
    % Map Human Proteins Given
    h_count=length(Vp_H_Seq);
    H_mapedSeq=zeros(h_count,length(AAs));
    for hin=1:h_count
        CrntSeq=Vp_H_Seq{hin,1}; % Character Array
        for n=1:length(AAs)
            H_mapedSeq(hin,n)=(length(strfind(CrntSeq,AAs(n))))/length(CrntSeq);
        end
    end
    % Assigning variable names to match previously built code
    F_HSeq=H_mapedSeq;
    F_HSeq_norm=F_HSeq;
    
    % Combine feature Vectors (Vp,Hp)
    % SVM data entry: row = sample, col= feature
    FeatureVectors=[repmat(F_Vseq_norm,size(F_HSeq,1),1),F_HSeq_norm];
end

%% Descriptor 3: Dipeptide Composition
if strcmp(Descriptor,'Dipeptide_Composition')
    %% Normalized Frequency of each Dipeptides in each protein
    AAs='GPAVLIMCFYWHKRQNEDST';
    % The single Virus Protein given
    Vp_mapedSeq=[];
    CrntVpSeq=V_Seq{1,Vp}; %Vp is a functin input of virus protein Index
    CrntVp_mapedSeq=zeros(1,length(AAs)^2);
    poin=0;
    for n=1:length(AAs)
        for n2=1:length(AAs)
            poin=poin+1;
            CrntDi=[AAs(n) AAs(n2)];
            CrntVp_mapedSeq(1,poin)=(length(strfind(CrntVpSeq,CrntDi)))/(length(CrntVpSeq)-1);
        end
    end
    % Assigning variable names to match previously built code
    Vp_mapedSeq=CrntVp_mapedSeq;
    F_Vseq=Vp_mapedSeq;
    F_Vseq_norm=F_Vseq;
    
    % Map Human Proteins Given
    h_count=length(Vp_H_Seq);
    H_mapedSeq=zeros(h_count,length(AAs)^2);
    for hin=1:h_count
        CrntSeq=Vp_H_Seq{hin,1}; % Character Array
        poin=0;
        for n=1:length(AAs)
            for n2=1:length(AAs)
                poin=poin+1;
                CrntDi=[AAs(n) AAs(n2)];
                H_mapedSeq(hin,n)=(length(strfind(CrntSeq,CrntDi)))/(length(CrntSeq)-1);
            end
        end
    end
    % Assigning variable names to match previously built code
    F_HSeq=H_mapedSeq;
    F_HSeq_norm=F_HSeq;
    
    % Combine feature Vectors (Vp,Hp)
    % SVM data entry: row = sample, col= feature
    FeatureVectors=[repmat(F_Vseq_norm,size(F_HSeq,1),1),F_HSeq_norm];
    
end

end

function Fnorm=NormalizeSVM(Norm_Type,F)
if strcmp(Norm_Type,'Cui')
    Fnorm=((F-min(F))/exp(max(F)-min(F)))-1;
elseif strcmp(Norm_Type,'Shen')
    Fnorm=(F-min(F))/(max(F)-min(F));
elseif strcmp(Norm_Type,'None')
    Fnorm=F;
end
end