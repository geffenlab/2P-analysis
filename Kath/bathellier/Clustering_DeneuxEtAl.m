
%% Clustering analysis
clear; close all;

% Loading the trial averaged responses of each neuron (DataM), there adress (IDX: which mouse and recording they were taken from), their signal to noise ratio, 
load('ClusteringDataNeuropilAllSounds.mat','DataM','IDX','SNR','mouse','Stims','ANOV')

% Parameters
dt=0.031731; % Frame interval in (s)
Nbefore=32;  % Number of frame before stimulus onset
t=([2:109]-Nbefore)*dt; % Time vector
Nc=50;   % Max number of cluster (dataset dependant, just try out different values)
SizeMin=10; % Minimum number of cells in a cluster (idem)
CorrThresh=0.1; % Threshold for re-assignement of low SNR cells to detected clusters


% Do some filtering first (Gaussian with one frame STD)
DataM=gaussianfilter(DataM,1,1,1);


% This is to remove cells which do not respond significantly to anything
idx=ANOV<0.05/length(SNR); NsignNeuron=size(idx)
DataM=DataM(idx,:,:); IDX=IDX(idx,:); StD=StD(idx); SNR=SNR(idx); ANOV=ANOV(idx);

% Size of remaining data
[m n p]=size(DataM);

% This is to select the 20% of cells with best SNR in each recorded population. You will need to adapt that to your data.
f=0.2; 
snr_top=[]; snr_bot=[];
for i=1:5
    for j=1:6
        idx=find(IDX(:,2)==i & IDX(:,3)==j);
        snr=SNR(idx);
        [tmp ix]=sort(snr,'descend');
        snr_top=[snr_top ; idx(ix(1:round(length(ix)*f)))];
        snr_bot=[snr_bot ; idx(ix(round(length(ix)*f)+1:end))];
    end
end
M=length(snr_top);
snr_idx=[snr_top; snr_bot];

% Compute correlation matrix
DataSel=reshape(DataM(snr_idx(1:M),30:end,:),[M (n-29)*p]);
DataSel=DataSel./repmat(std(DataSel,1,2),[1 (n-29)*p]);

% Compute correlation matrix
CC=pdist(DataSel); 

% Hierarchical clustering
Z=linkage(CC,'complete'); 

% Organize and plot the clustering
figure;
[H,T,perm] = dendrogram(Z,0); 

% Select a maximum of 50 cluster (this number is dataset dependant)
H=cluster(Z,'maxClust',Nc);
Hp=H(perm);

% Calculate the size (ncell) and mean response (Data) of each cluster 
idx=snr_idx(1:M);
for i=1:Nc
    is(i)=find(Hp==i,1,'first'); 
    ie(i)=find(Hp==i,1,'last');
    ncell(i)=ie(i)-is(i)+1;
    Data(i,:,:)=mean(DataM(idx(H==i),:,:),1);
end

% remove clusters that are too small
Data=Data(ncell>SizeMin,:,:);

% Then you can plot that
PlotMultipleTraces(permute(Data,[1 3 2]),0.0315,0.5,0.02);

% Here I recommand to look at all cluster and remove those that are obviously noise. Then you can proceed to the next step (this is dataset dependant so you have to implement that).

% Agglomeration based on correlation coeff. for cluster
IdxOut=[idx(Iout); snr_idx(M+1:end)]; 
cc=corrcoef(reshape([DataM(IdxOut,30:end,:) ; Data(:,30:end,:)],[numel(IdxOut)+nc (n-29)*p])');
cc=cc(end-nc+1:end,1:end-nc);
[tmp, Cidx]=max(cc,[],1);
Cidx(tmp<CorrThresh)=0;


for i=1:nc
    IdxIn2{i}=[IdxIn{i}; IdxOut(Cidx==i)]';
    IdxIn{i}=IdxIn{i}';
    Ncells0(i)=numel(IdxIn{i});
    Data(i,:,:)=mean(DataM(IdxIn2{i},:,:),1);
    Ncells(i)=numel(IdxIn2{i});
end
IdxOut=IdxOut(Cidx==0)';


% Then you can plot that again
PlotMultipleTraces(permute(Data,[1 3 2]),0.0315,0.5,0.02);




