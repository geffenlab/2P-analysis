%% Make averages
mouse(1).base='Louise';
mouse(1).dates={...
'140417', ...
'140418', ...
'140429am', ...
'140429pm',...
'140507'};

slist_l=[1:6 11 12:17 22 25:32 35:42 45:52 55:62];
 slist_l=[12:17 22 35:42 55:62];  % White noise only

mouse(2).base='jason';
mouse(2).dates={...
'150327',...
'150328',...
'150330',...
'150401'};

mouse(3).base='kelly';
mouse(3).dates={...
'150327',...
'150330',...
'150331',...
'150403'};

slist_jk=[1:14 17:24 27:34 37:44 47:54];
 slist_jk=[8:14 27:34 47:54]; % White noise only

mouse(4).base='m4';
mouse(4).dates={...
'141128'};

mouse(5).base='Dot0';
mouse(5).dates={...
    '150324'};

slist_md=[41:45 47:48 1:2:8 13:2:20 2:2:8 14:2:20]; % White noise only
 
DT=0.031731;
Nbefore=32;
t=([1:109]-Nbefore)*DT;

Tlist=[3 7 10 14 22 38 30 46];
tau=2;

M=0; SNR=[]; DataM=[]; IDX=[]; ANOV=[]; StD=[]; T=[]; n4=20; TTES=[];
for i=1:3%length(mouse)
    for j=1:length(mouse(i).dates)
        file=['C:\Data Drive\Calcium_Imaging\Data\2016_DataRamp-Neuropile\' mouse(i).dates{j} '_' mouse(i).base '\data.mat']
        load(file) 
        m=size(Traces,1);
        IDXtmp=[(1:m)' i*ones(m,1) j*ones(m,1)];
        if(i==1)
            idx=~isnan(mean(Traces(:,:),2)); 
            Traces=Traces(idx,:,:,:); IDXtmp=IDXtmp(idx,:);
            [m, n, p, q]=size(Traces);
            TMP=zeros([m, 109, p, q]);
            tl=([1:n]-Nbefore)*dt;
            Traces=permute(interp1(tl',permute(Traces,[2 1 3 4]),t'),[2 1 3 4]);
            Traces=Traces(:,:,slist_l,:);
            Stims=Stims(slist_l);
        elseif(i<4)
            idx=~isnan(mean(Traces(:,:),2));
            Traces=Traces(idx,:,:,:); IDXtmp=IDXtmp(idx,:);
            Traces=Traces(:,:,slist_jk,:);
            Stims=Stims(slist_jk);
        else
            idx=~isnan(mean(Traces(:,:),2));
            Traces=Traces(idx,:,:,:); IDXtmp=IDXtmp(idx,:);
            Traces=Traces(:,:,slist_md,:);
            Stims=Stims(slist_md);
        end
        Traces=Traces(:,2:end,:,:);
        [m, n, p, q]=size(Traces);
        M=M+m;
        Traces=[diff(Traces,1,2) zeros(m,1,p,q)]+[Traces(:,1:end-1,:,:)  zeros(m,1,p,q)]*DT/tau;
        Traces=Traces-repmat(mean(mean(Traces(:,1:31,:,:),2),4),[1 n 1 q]);
        SNR=[SNR ; mean(mean(abs(mean(Traces,4)),2),3)./mean(mean(std(Traces,1,4),2),3)];
        panov=zeros(1,m); test=zeros(1,m);
        for k=1:m
            panov(k) = kruskalwallis(squeeze(mean(Traces(k,31:end,:,:),2))',[],'off');
            [~, test(k)]=signrank(squeeze(mean(Traces(k,1:31,:),2))-squeeze(mean(Traces(k,32:end,:),2)));
        end
        ANOV=[ANOV; panov']; TTES=[TTES; test'];
        DataM=[DataM; mean(Traces,4)];
        StD=[StD; std(Traces(:,:),1,2)];
        IDX=[IDX; IDXtmp];
%         if(isempty(T))
%             T=Traces(:,:,Tlist,:);
%         else
%             n4=size(T,4);
%             T=[T(:,:,:,1:min(n4,size(Traces,4)));Traces(:,:,Tlist,1:min(n4,size(Traces,4)))];
%         end
    end
end


%save('Traces6085CstAndRamp2s.mat','T','StD')
save('ClusteringDataNeuropilWN.mat','DataM','IDX','SNR','mouse','Stims','StD','ANOV','TTES')

%% Load averages
clear; close all;
load('ClusteringDataNeuropilWN.mat','DataM','IDX','SNR','mouse','Stims','StD','ANOV')
dt=0.031731; Nbefore=32; t=([2:109]-Nbefore)*dt; offset=0.001;

idx=StD<1;
DataM=DataM(idx,:,:); IDX=IDX(idx,:); StD=StD(idx,:); SNR=SNR(idx,:); ANOV=ANOV(idx,:);
DataMf=gaussianfilter(DataM,1,1.5,1);
f=0;
for i=1:5
    for j=1:5
        idx=(IDX(:,2)==i & IDX(:,3)==j);
        if(sum(idx))
            f=f+1;
            DataMfG(f,:,:)=mean(DataMf(idx,:,:),1);
        end
    end
end

dat=mean(DataMf(:,:,[3 7 13 15 21 23]),1);
save('UpDown.mat','dat');

idx=t>-0.25;
figure(1);
subplot(5,1,1)
fn_errorbar(t(idx),squeeze(mean(DataMf(:,idx,[15 23]),1))+offset,squeeze(std(DataMfG(:,idx,[15 23]),1,1))/sqrt(14),'patch')
subplot(5,1,2)
fn_errorbar(t(idx),squeeze(mean(DataMf(:,idx,[9:2:15]),1))+offset,squeeze(std(DataMfG(:,idx,[9:2:15]),1,1))/sqrt(14),'patch')
subplot(5,1,3)
fn_errorbar(t(idx),squeeze(mean(DataMf(:,idx,[17:2:23]),1))+offset,squeeze(std(DataMfG(:,idx,[17:2:23]),1,1))/sqrt(14),'patch')

% Integrale differences
Iup=[9 11 13 15]; Idown=[17 19 21 23];
Toff=[0.1 0.25 1 2];

for i=1:4
    IntUp60(:,i)=sum(DataMfG(:,t>=0&t<=Toff(i)+0.45,Iup(i))+offset,2);
    IntDo60(:,i)=sum(DataMfG(:,t>=0&t<=Toff(i)+0.45,Idown(i))+offset,2);
    IntUp50(:,i)=sum(DataMfG(:,t>=0&t<=Toff(i)+0.45,Iup(i)-1)+offset,2);
    IntDo50(:,i)=sum(DataMfG(:,t>=0&t<=Toff(i)+0.45,Idown(i)-1)+offset,2);
end

figure;
subplot(2,1,1)
xx=(IntUp60-IntDo60);
yy=(IntUp60+IntDo60)*0.5;
errorbar(Toff,mean(xx,1)./mean(yy,1),std(xx,1,1)./mean(yy,1)/sqrt(14))
for i=1:4
    [~, pval(1,i)]=kstest(IntUp60(:,i)-IntDo60(:,i));
    pvals(1,i)=signrank(IntUp60(:,i)-IntDo60(:,i),0,'tail','right');
end
subplot(2,1,2)
xx=(IntUp50-IntDo50);
yy=(IntUp50+IntDo50)*0.5;
errorbar(Toff,mean(xx,1)./mean(yy,1),std(xx,1,1)./mean(yy,1)/sqrt(14))
for i=1:4
    [~, pval(2,i)]=kstest((IntUp50(:,i)-IntDo50(:,i)));
    pvals(2,i)=signrank(IntUp50(:,i)-IntDo50(:,i),0,'tail','right');
end



load('ClusteringDataNeuropilAllSounds.mat','DataM','IDX','SNR','mouse','Stims','StD','ANOV')
dt=0.031731; Nbefore=32; t=([2:109]-Nbefore)*dt;

idx=StD<1;
DataM=DataM(idx,:,:); IDX=IDX(idx,:); StD=StD(idx,:); SNR=SNR(idx,:); ANOV=ANOV(idx,:);
DataMf2=gaussianfilter(DataM,1,1.5,1);


f=0;
for i=1:3
    for j=1:5
        idx=(IDX(:,2)==i & IDX(:,3)==j);
        if(sum(idx))
            f=f+1;
            DataMfG2(f,:,:)=mean(DataMf2(idx,:,:),1);
        end
    end
end

figure(1);
idx=t>-0.25; idx=idx(1:end-3);
subplot(5,1,4)
fn_errorbar(t(idx),squeeze(mean(DataMf2(:,idx,[22 38]),1))+offset,squeeze(std(DataMfG2(:,idx,[22 38]),1,1))/sqrt(14),'patch')

figure;
subplot(2,1,1)
fn_errorbar(t(idx),squeeze(mean(DataMf2(:,idx,[16 18 20 22]),1))+offset,squeeze(std(DataMfG2(:,idx,[16 18 20 22]),1,1))/sqrt(14),'patch')
subplot(2,1,2)
fn_errorbar(t(idx),squeeze(mean(DataMf2(:,idx,[32 34 36 38]),1))+offset,squeeze(std(DataMfG2(:,idx,[32 34 36 38]),1,1))/sqrt(14),'patch')


% Integrale differences
Iup=[16 18 20 22]; Idown=[32 34 36 38];
Toff=[0.1 0.25 1 2];
for i=1:4
    IntUp60c(:,i)=sum(DataMfG2(:,t>=0&t<=Toff(i)+0.45,Iup(i))+offset,2);
    IntDo60c(:,i)=sum(DataMfG2(:,t>=0&t<=Toff(i)+0.45,Idown(i))+offset,2);
    IntUp50c(:,i)=sum(DataMfG2(:,t>=0&t<=Toff(i)+0.45,Iup(i)-1)+offset,2);
    IntDo50c(:,i)=sum(DataMfG2(:,t>=0&t<=Toff(i)+0.45,Idown(i)-1)+offset,2);
end

figure;
subplot(2,1,1)
hold on
xx=(IntUp60c-IntDo60c);
yy=(IntUp60c+IntDo60c)*0.5;
errorbar(Toff,mean(xx,1)./mean(yy,1),std(xx,1,1)./mean(yy,1)/sqrt(12))
hold off
for i=1:4
    [~, pval(3,i)]=kstest(IntUp60c(:,i)-IntDo60c(:,i));
    pvals(3,i)=signrank(IntUp60c(:,i)-IntDo60c(:,i),0,'tail','right');
end
subplot(2,1,2)
xx=(IntUp50c-IntDo50c);
yy=(IntUp50c+IntDo50c)*0.5;
hold on
errorbar(Toff,mean(xx,1)./mean(yy,1),std(xx,1,1)./mean(yy,1)/sqrt(14))
hold off
for i=1:4
    [~, pval(4,i)]=kstest(IntUp50c(:,i)-IntDo50c(:,i));
    pvals(4,i)=signrank(IntUp50c(:,i)-IntDo50c(:,i),0,'tail','right');
end

% Intensity responses
IdxWN=[8:14]; IdxCpx=[1:7];
dB=[50 55 60 65 70 80 85];
Toff=0.25;

for i=1:7
   IntWN(:,i)=sum(DataMfG2(:,t>=0&t<=Toff+0.45,IdxWN(i))+offset,2);
   IntCpx(:,i)=sum(DataMfG2(:,t>=0&t<=Toff+0.45,IdxCpx(i))+offset,2);
end

figure;
subplot(2,1,1)
hold on


errorbar([dB;dB]',[mean(IntWN,1); mean(IntCpx,1)]',[std(IntWN,1,1); std(IntCpx,1,1)]'/sqrt(12))


%% Clustering analysis
clear; close all;
load('ClusteringDataNeuropilAllSounds.mat','DataM','IDX','SNR','mouse','Stims','StD','ANOV')
dt=0.031731; Nbefore=32; t=([2:109]-Nbefore)*dt;

idx=StD<1; Nneuron=size(idx)
DataM=DataM(idx,:,:); IDX=IDX(idx,:); StD=StD(idx,:); SNR=SNR(idx,:); ANOV=ANOV(idx,:);
DataMf=gaussianfilter(DataM,1,1.5,1);


% Compute correlation matrix 
DataM=gaussianfilter(DataM,1,1,1);
f=0.2; 


idx=ANOV<0.05/length(SNR); NsignNeuron=size(idx)
DataM=DataM(idx,:,:); IDX=IDX(idx,:); StD=StD(idx); SNR=SNR(idx); ANOV=ANOV(idx);

[m n p]=size(DataM);

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

% M=round(f*m);
% [tmp snr_idx]=sort(SNR,'descend');
DataSel=reshape(DataM(snr_idx(1:M),30:end,:),[M (n-29)*p]);
% CC=corrcoef(DataSel');
% [T,perm] = dendrogramNOPLOT(linkage(squareform(1-CC-eye(M).*(1-CC)),'complete'),0); % complete, average, single(deflt), weighted
DataSel=DataSel./repmat(std(DataSel,1,2),[1 (n-29)*p]);
CC=pdist(DataSel);
Z=linkage(CC,'complete');  % complete, average, single(deflt), weighted
% [T,perm] = dendrogramNOPLOT(Z,0); % complete, average, single(deflt), weighted
figure;
[H,T,perm] = dendrogram(Z,0);

Nc=50;
SizeMin=10;
CorrThresh=0.1;
H=cluster(Z,'maxClust',Nc);
Hp=H(perm);

idx=snr_idx(1:M);
for i=1:Nc
    is(i)=find(Hp==i,1,'first'); 
    ie(i)=find(Hp==i,1,'last');
    ncell(i)=ie(i)-is(i)+1;
    Data(i,:,:)=mean(DataM(idx(H==i),:,:),1);
end

rc=[1:7  15:22 31:38];
rw=[8:14 23:30 39:46];
PlotMultipleTraces(permute(Data(ncell>SizeMin,:,rc),[1 3 2]),0.0315,0.5,0.04);
PlotMultipleTraces(permute(Data(ncell>SizeMin,:,rw),[1 3 2]),0.0315,0.5,0.04);

CC=squareform(CC);
CCclust=CC(:,perm);
CCclust=CCclust(perm,:);
figure; imagesc(CCclust); title(['> SizeMin = ' num2str(sum(ncell>SizeMin))])
nc=0; Iout=zeros(size(H));
for i=1:Nc
   line([is(i) ie(i) ie(i) is(i) is(i)],[is(i) is(i) ie(i) ie(i) is(i)],'Color','k')
   if(ncell(i)>SizeMin)
       nc=nc+1;
       text(ie(i),is(i),num2str(nc))  
       IdxIn{nc}=idx(H==i);
   else
       Iout=Iout|(H==i);
   end
end

% Organization and selection of clusters 
permc=[[4 5 6 7 8 14 17 18 20 21 23 24 25] 1 2 3 9 10 11 12 13 15 16 19 22 26];
figperm=[[9 10 6 13 12 11 5 7 4 2 8 1 3] 14:26]; 
permc=permc(figperm);



% Agglomeration based on correlation coeff.
Data=Data(ncell>SizeMin,:,:);
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



PlotMultipleTraces(permute(Data(:,:,rc),[1 3 2]),0.0315,0.5,0.02);
PlotMultipleTraces(permute(Data(:,:,rw),[1 3 2]),0.0315,0.5,0.02);



%%%
nr=13;
Dataf=Data(permc(1:end-nr),:,:);
Ncellsf=Ncells(permc(1:end-nr));
Ncells0=Ncells0(permc(1:end-nr));
IdxInf=IdxIn2(permc(1:end-nr));
IdxIn=IdxIn(permc(1:end-nr));
IdxOut=[IdxOut [IdxIn2{permc(end-nr+1:end)}]];

% plot selected clusters in PCA space
dat=reshape(DataM([IdxIn{:}],30:end,:),sum(Ncells0),(n-29)*p);
dat=dat./repmat(std(dat,1,2),1,(n-29)*p);
[tmp, PCs, spect]=princomp(dat);

colors = [
    150 50 255
    233 217 143
    255 168 211
    250 120 0
    0   0   0
    0   0   255
    0   100 160
    0   255 0
    0   138 0
    191 0   191
    255 255 0
    200 0 0
    255 0   0
    ]/255;

figure; 
Nce=[1 cumsum(Ncells0)];
for i=1:length(Ncells0)-1
    hold on
    plot3(PCs(Nce(i):Nce(i+1),1),PCs(Nce(i):Nce(i+1),2),PCs(Nce(i):Nce(i+1),3),'.','Color',colors(i,:));
    hold off
end


% % Cluster removal
% for i=[1]
%    Dataf(i,:,:)=(Dataf(i,:,:)*Ncellsf(i)+Dataf(i+1,:,:)*Ncellsf(i+1))/(sum(Ncellsf(i:i+1)));
%    Dataf=Dataf([1:i i+2:end],:,:);
%    Ncellsf(i)=Ncellsf(i)+Ncellsf(i+1);
%    Ncellsf=Ncellsf([1:i i+2:end]);
%    IdxInf{i}=[IdxInf{i} IdxInf{i+1}];
%    IdxInf=IdxInf([1:i i+2:end]);
% end

Ncellsf
PlotMultipleTraces(permute(Dataf(:,20:end,rc),[1 3 2]),0.0315,0.5,0.04);
PlotMultipleTraces(permute(Dataf(:,20:end,rw),[1 3 2]),0.0315,0.5,0.04);

dcomplex=Dataf(:,:,rc);
whitenoise=Dataf(:,:,rw);
weigths=[Ncellsf numel(IdxOut)];


%%% Correlation matrix: all cells
sels=1:46; %[2 7 10 14 43:46 27:30 35:38 19:22];
NselC=length(permc)-nr; p1=length(sels);

DataMff=reshape(DataM([IdxInf{:}],30:end,sels),[numel([IdxInf{:}]) (n-29)*p1])';
ccf=corrcoef(DataMff);

% DataMff=DataMff'./repmat(std(DataMff',1,2),[1 (n-29)*p]);
% ccf=squareform(pdist(DataMff));

mcc=mean(ccf,1);



label=zeros(size(IDX,1),1);

permk=[]; norm=[];
for i=1:NselC
    is(i)=length([IdxInf{1:i-1}])+1;
    ie(i)=length([IdxInf{1:i}]);
    list=is(i):ie(i);  
    [tmp idx]=sort(mcc(list),'ascend'); 
    permk=[permk list(idx)];
    tmp=mean(ccf(list,list),1);
    norm=[norm 1-tmp];
end
ccf=ccf(permk,:); ccf=ccf(:,permk);
tmp=repmat(norm,[length(norm) 1]);
norm=(tmp+tmp')/2;
norm=norm(permk,:); norm=norm(:,permk);

% figure; imagesc(1-ccf)
figure; imagesc((1-ccf)./norm)


for i=1:NselC
    line([is(i) ie(i) ie(i) is(i) is(i)],[is(i) is(i) ie(i) ie(i) is(i)],'Color','k')
    % text(ie(i),is(i),num2str(i))
    label(IdxInf{i})=i;
end
save('C:\Data Drive\Calcium_Imaging\Data\briceclustersNeuropil.mat','dcomplex','whitenoise','weigths','IdxInf','IDX','label','DataM');


% %%% Correlation matrix: 20% SNR
% ccf=corrcoef(reshape(DataM([IdxIn{:}],30:end,:),[numel([IdxIn{:}]) (n-29)*p])');
% figure; imagesc(1-ccf)
% NselC=length(permc)-4;
% for i=1:NselC
%    is(i)=length([IdxIn{1:i-1}])+1;
%    ie(i)=length([IdxIn{1:i}]);
%    line([is(i) ie(i) ie(i) is(i) is(i)],[is(i) is(i) ie(i) ie(i) is(i)],'Color','k')
% end


save('DatForThomasNeuropil.mat','mouse','IDX','label','ccf','IdxInf')


%%% PCA plots
ru=30; rd=46;
[tmp, PCs]=princomp(reshape(DataM(:,:,[ru rd]),[size(DataM,1) size(DataM,2)*2])');
dec=10;
PCs=reshape(PCs',[size(DataM,1) size(DataM,2) 2]);
PCs=permute(interp1(1:size(PCs,2),permute(PCs,[2 1 3]),1:1/dec:size(PCs,2),'spline'),[2 1 3]);
figure; plot3(PCs(1,:,1),PCs(2,:,1),PCs(3,:,1),'r',PCs(1,:,2),PCs(2,:,2),PCs(3,:,2),'b',PCs(1,1,1),PCs(2,1,1),PCs(3,1,1),'.r',PCs(1,1,2),PCs(2,1,2),PCs(3,1,2),'.b')

%%%% Plot integrals
Tin=[0.25*ones(1,7) 0.1 0.1 0.25 0.25 1 1 2 2];
Tin=[Tin(1:7) Tin(8:end) Tin(8:end)];

% find populations
idx=[IdxInf{:}];
popcode=IDX(:,2)+10*IDX(:,3);
pop=unique(popcode);
pc=popcode(idx);
    
for i=1:length(Tin)
    for j=1:length(pop)
        IntD(j,i)=sum(sum(DataM(pc==pop(j),32+(1:round(Tin(i)/dt)+13),rw(i))));
    end
end

figure;
subplot(2,1,1)
errorbar([0.1 0.25 1 2],2*mean(IntD(:,8:2:15)-IntD(:,16:2:23),1)./mean(IntD(:,8:2:15)+IntD(:,16:2:23),1),2*std(IntD(:,8:2:15)-IntD(:,16:2:23),1,1)./mean(IntD(:,8:2:15)+IntD(:,16:2:23),1)/sqrt(12))
xlim([0.05 2.5])
subplot(2,1,2)
errorbar([0.1 0.25 1 2],2*mean(IntD(:,9:2:15)-IntD(:,17:2:23),1)./mean(IntD(:,9:2:15)+IntD(:,17:2:23),1),2*std(IntD(:,9:2:15)-IntD(:,17:2:23),1,1)./mean(IntD(:,9:2:15)+IntD(:,17:2:23),1)/sqrt(12))
xlim([0.05 2.5])

% Integral diff for each cluster
for i=1:length(Tin)
        IntCW(:,i)=sum(Dataf(:,32+(1:round(Tin(i)/dt)+13),rw(i)),2);
end
for i=1:length(Tin)
        IntCC(:,i)=sum(Dataf(:,32+(1:round(Tin(i)/dt)+13),rc(i)),2);
end

figure; 
subplot(2,2,1)
imagesc((IntCW(:,9:2:15)-IntCW(:,17:2:end))',[-.7 , 0.7]); title('WN 60-85')
subplot(2,2,2)
imagesc((IntCW(:,8:2:15)-IntCW(:,16:2:end))',[-.7 , 0.7]); title('WN 50-85')
subplot(2,2,3)
imagesc((IntCC(:,9:2:15)-IntCC(:,17:2:end))',[-.7 , 0.7]); title('8khz 60-85')
subplot(2,2,4)
imagesc((IntCC(:,8:2:15)-IntCC(:,16:2:end))',[-.7 , 0.7]); title('8khz 50-85')

colormap(spring)

%% Load averages 2
clear; 
% close all;

rc=[1:7  15:22 31:38];
rw=[8:14 23:30 39:46];
f=0.2; 

load('ClusteringDataWN.mat','DataM','IDX','SNR','mouse')
DataM1=DataM; SNR1=SNR; IDX1=IDX;

load('ClusteringData.mat','DataM','IDX','SNR','mouse')
DataM=[DataM(:,:,rw);DataM1]; SNR=[SNR; SNR1]; IDX=[IDX; IDX1];
[m n p]=size(DataM);

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

% M=round(f*m);
% [tmp snr_idx]=sort(SNR,'descend');

% Compute correlation matrix 
DataM=gaussianfilter(DataM,1,1,1);  %  1


DataSel=reshape(DataM(snr_idx(1:M),30:end,:),[M (n-29)*p]);
% CC=corrcoef(DataSel');
% [T,perm] = dendrogramNOPLOT(linkage(squareform(1-CC-eye(M).*(1-CC)),'complete'),0); % complete, average, single(deflt), weighted
DataSel=DataSel./repmat(std(DataSel,1,2),[1 (n-29)*p]);
CC=pdist(DataSel);
Z=linkage(CC,'complete');  % complete  %%% average, single(deflt), weighted
% [T,perm] = dendrogramNOPLOT(Z,0); % complete, average, single(deflt), weighted
[H,T,perm] = dendrogram(Z,0); 

Nc=50;            % 50
SizeMin=10;       % 10
CorrThresh=0.1;   % 0.1 
H=cluster(Z,'maxClust',Nc);
Hp=H(perm);

idx=snr_idx(1:M);
for i=1:Nc
    is(i)=find(Hp==i,1,'first'); 
    ie(i)=find(Hp==i,1,'last');
    ncell(i)=ie(i)-is(i)+1;
    Data(i,:,:)=mean(DataM(idx(H==i),:,:),1);
end


PlotMultipleTraces(permute(Data(ncell>SizeMin,:,:),[1 3 2]),0.0315,0.5,0.04);

CC=squareform(CC);
CCclust=CC(:,perm);
CCclust=CCclust(perm,:);
figure; imagesc(CCclust); title(['> SizeMin = ' num2str(sum(ncell>SizeMin))])
nc=0; Iout=zeros(size(H));
for i=1:Nc
   line([is(i) ie(i) ie(i) is(i) is(i)],[is(i) is(i) ie(i) ie(i) is(i)],'Color','k')
   if(ncell(i)>SizeMin)
       nc=nc+1;
       text(ie(i),is(i),num2str(nc))  
       IdxIn{nc}=idx(H==i);
   else
       Iout=Iout|(H==i);
   end
end

% Agglomeration based on correlation coeff.
Data=Data(ncell>SizeMin,:,:);
IdxOut=[idx(Iout); snr_idx(M+1:end)]; 
cc=corrcoef(reshape([DataM(IdxOut,30:end,:) ; Data(:,30:end,:)],[numel(IdxOut)+nc (n-29)*p])');
cc=cc(end-nc+1:end,1:end-nc);
[tmp, Cidx]=max(cc,[],1);
Cidx(tmp<CorrThresh)=0;

for i=1:nc
    IdxIn{i}=[IdxIn{i}; IdxOut(Cidx==i)];
    Data(i,:,:)=mean(DataM(IdxIn{i},:,:),1);
    Ncells(i)=numel(IdxIn{i});
end
IdxOut=IdxOut(Cidx==0);

PlotMultipleTraces(permute(Data,[1 3 2]),0.0315,0.5,0.02);


% %%% Cluster reduction
permc=[1 2 [3 12 13] 8 5 18];
Dataf=Data(permc(1:8),:,:);
Ncellsf=Ncells(permc(1:8));
IdxInf=IdxIn(permc(1:end-1));
IdxOut=[IdxOut ; IdxIn{permc(end)}];

% for i=[4 3]
%    Dataf(i,:,:)=(Dataf(i,:,:)*Ncellsf(i)+Dataf(i+1,:,:)*Ncellsf(i+1))/(sum(Ncellsf(i:i+1)));
%    Dataf=Dataf([1:i i+2:end],:,:);
%    Ncellsf(i)=Ncellsf(i)+Ncellsf(i+1);
%    Ncellsf=Ncellsf([1:i i+2:end]);
% end

Ncellsf
PlotMultipleTraces(permute(Dataf(:,:,:),[1 3 2]),0.0315,0.5,0.02);

% %% 
% clear all; 
% load('Traces6085CstAndRamp2s.mat','T','StD')
% idx=StD<1;
% T=T(idx,:,:,:);
% dt =0.03173; 
% list.baseline=32*dt;
% names={'Cplx60dB' 'Cplx85dB' 'WN60dB' 'WN85dB' 'CplxUp' 'CplxDown' 'WNUp' 'WNDown' };
% sounds=[3 4 7 8];
% TimeCorrelation2(T,names,dt,sounds,20:108);
% sounds=[3 4 7 8]-2;
% TimeCorrelation2(T,names,dt,sounds,20:108);
% sounds=5:8;
% TimeCorrelation2(T,names,dt,sounds,20:108);
% 
