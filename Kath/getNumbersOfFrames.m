function n = getNumbersOfFrames(folderNames)

for ii=1:length(folderNames)
    ts = dir([folderNames{ii} '*.tif']);
    for jj=1:length(ts)      
        fname = [folderNames{ii} ts(jj).name];
        ts_info = imfinfo(fname);
        num_images(jj) = numel(ts_info);
    end
    n(ii) = sum(num_images);
    num_images=[];
end








