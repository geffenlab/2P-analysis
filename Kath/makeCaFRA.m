function fra = makeCaFRA(raster,window,tones,attns,method)

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
        
        switch method
            case 'trapz'
                auc(ii,:) = trapz(psth(window(1):window(2),:));       
            case 'peak'
                auc(ii,:) = max(psth(window(1):window(2),:));    
        end
    end
    fra(:,:,jj) = auc';
end