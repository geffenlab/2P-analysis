% z-motion analysis


% load the Z-stack
path_to_Zfile = 'C:\data\K024\20161026K024_tifStacks\8\ZSeries-10262016-30micron162.5-016_1.tif';
sframe = 1;
num2read = 10000;
zstack=bigread2_KCW(path_to_Zfile,sframe,num2read);
zstack = zstack(:,:,11:end);



% load the time stack
path_to_file = 'C:\data\K024\20161026K024_tifStacks\3\TSeries-10262016-rang-regStim-011_2_Z.tif';
sframe = 1;
num2read = 10000;
imData=bigread2_KCW(path_to_file,sframe,num2read);

% Correlate each plane in time series with the z stack
for jj = 1:size(zstack,3)
    parfor ii = 1:size(imData,3)
        r = corr2(imData(206:306,206:306,ii),zstack(206:306,206:306,jj));
        R(ii,jj) = r;
    end
end