function done = tiffStackMaker(folder,identifier,saveFileName,saveFolder)

fnames = dir([folder '*' identifier '*.tif']);
if isempty(fnames)
    fnames = dir([folder '*' identifier '*.tiff']);
    if isempty(fnames)
        fnames = dir([folder '*' identifier '.tif']);
        if isempty(fnames)
            folder = [folder '\'];
            fnames = dir([folder '*' identifier '.tif']);
            %                 folder = folder(1:end-1);
        else
            disp('Can''t find tiffs!')
            return
        end
    end
end


% stackInfo = imfinfo([folder fnames(1).name]);
if ~isdir(saveFolder)
    mkdir(saveFolder)
end
% stack = imread([folder fnames(1).name]);
ind = 1; nFiles = 0;
% imwrite(stack, [saveFolder saveFileName '_' num2str(ind) '.tif'], 'writemode', 'append','Compression','None');
for k = 1:length(fnames)
    nFiles = nFiles+1;
    stack = imread([folder fnames(k).name ]);
%     if ~exist([saveFolder saveFileName '_' sprintf('%02d',ind) '.tif'],'file')
        imwrite(stack, [saveFolder saveFileName '_' sprintf('%02d',ind) '.tif'], 'writemode', 'append','Compression','None');
%     end
        if mod(k,1000)==0
            disp(['Written ' num2str(k) ' frames of ' num2str(length(fnames))])
        end
        if nFiles==2000
            ind = ind+1; disp(ind)
            nFiles = 0;
        end
end
disp('Written all tifs')
done = 1;


