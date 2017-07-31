% Add Bathellier SNR
clear
dataLoc = 'E:\dataAnalysed\Kath\Bathellier\';
files = dir([dataLoc '*.mat']);

load([dataLoc files(1).name])
raster = data.raster_DFF0;

% Add SNR
respWindow = data.nFramesPreStim+1:round(2*data.frameRate);

for ii=1:size(raster,3)
    
    SNR(ii) = trapz(abs(mean(raster(:,respWindow,ii),1))) /  sqrt(trapz(mean((raster(:,respWindow,ii) - mean(raster(:,respWindow,ii),1)).^2,1)))  ;
    
end

data.SNR = SNR;

save([dataLoc files(1).name],'data');