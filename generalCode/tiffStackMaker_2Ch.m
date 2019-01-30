function done = tiffStackMaker_2Ch(folder,identifier,saveFileName,saveFolder)

% Look for Tiffs
fnames = dir([folder '*' identifier '*.tif']);

if isempty(fnames)
    fnames = dir([folder '*' identifier '*.tiff']);
    if isempty(fnames)
        fnames = dir([folder '*' identifier '.tif']);
        if isempty(fnames)
            folder = [folder '\'];
            fnames = dir([folder '*' identifier '.tif']);
        else
            disp('Can''t find tiffs!')
            return
        end
    end
end


% addition for when there are two channels : sort by frame number rather
% than by channel to have channels interleaved in stack
a={fnames.name};
for ii=1:length(a)
s=strfind(a{ii},'_');
p=strfind(a{ii},'.ome.tif');
b(ii,:)=str2num(a{ii}(s(end)+1:p-1)); % find frame number
end
[B,I]=sort(b,1); % sorting by frame number
fnames_2ChOrd=fnames(I); % fnames ordered with Ch1 and Ch2 interleaved


% Get the scaling Factor (data is saved as 13 bits but matlab only has 16
% bitdepth therefore need to get the max across all the tifs and /max *
% 2^13 - I don't think we need to do this... the bit depth is 13 but the
% prairie view saves them as 16 so I think it is ok....?

if ~isdir(saveFolder)
    mkdir(saveFolder)
end


ind = 1; nFiles = 0;
for k = 1:length(fnames_2ChOrd)
    nFiles = nFiles+1;
    stack = imread([folder fnames_2ChOrd(k).name ]);
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


