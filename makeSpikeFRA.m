function fra = makeSpikeFRA(raster,window,tones,attns)

uT = unique(tones);
uA = unique(attns);

fra = zeros(length(uA),length(uT),size(raster,3));
for jj=1:size(raster,3)
    auc = zeros(length(uA),length(uT))';
    for ii=1:length(uT)
        psth=zeros(size(raster,2),length(uA));
        for kk=1:length(uA)
            rows = tones==uT(ii) & attns==uA(kk);
            psth(:,kk) = mean(raster(rows,:,jj));
        end
        
        auc(ii,:) = mean(psth(window(1):window(2),:));       
    end
    fra(:,:,jj) = auc';
end