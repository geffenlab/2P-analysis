function raster = makeCaRaster(traces,eventsOn,preOnsetTime,postOnsetTime,doDFoverF0) 
% traces (m x n matrix, where row is cell and column is time), eventsOn in same timescale as traces,preStimOnset in same timescale, do DFoverF0 is logical, if you want to subtract the mean and divide by the standard deviation of the preOnsetTime period


raster = zeros(length(eventsOn),postOnsetTime+preOnsetTime+2,size(traces,1));
for ii=1:length(eventsOn) % ignore 1st event because will be the start recording event - events from then on are relevant
    raster(ii,:,:) = traces(:,eventsOn(ii)-preOnsetTime-1: eventsOn(ii)+postOnsetTime)';
    if doDFoverF0==1
         raster(ii,:,:) = (squeeze(raster(ii,:,:))-mean(squeeze(raster(ii,1:preOnsetTime,:)),1))./std(squeeze(raster(ii,1:preOnsetTime,:)),[],1);
    end
end

