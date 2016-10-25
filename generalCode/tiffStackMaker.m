function done = tiffStackMaker(folder,identifier,saveFileName)

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


stackInfo = imfinfo([folder fnames(1).name]);

stack = imread([folder fnames(1).name]);
ind = 1; nFiles = 1;
imwrite(stack, [folder saveFileName '_' num2str(ind) '.tif'], 'writemode', 'append','Compression','None');
for k = 2:length(fnames)
    stack = imread([folder fnames(k).name ]);
    imwrite(stack, [folder saveFileName '_' num2str(ind) '.tif'], 'writemode', 'append','Compression','None');
    if mod(k,1000)==0
        disp(['Written ' num2str(k) ' frames of ' num2str(length(fnames))])
        toc
        tic
    end
    nFiles = nFiles+1;
    if nFiles+1>2000
        ind = ind+1; disp(ind)
        nFiles =1;
    end
end
disp('Written all tifs')
done = 1;


