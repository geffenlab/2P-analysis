function recInfo = readXmlFile(filename)
% read xml file
fclose('all');
x = fopen(filename);
sForm = '%s'; nSamples = 1;
stuffToExtract = {'.ome.tif','bitDepth','dwellTime','laserPower','laserWavelength','micronsPerPixel','objectiveLens','opticalZoom','pmtGain','ZAxis','Sequence type'};
ste = stuffToExtract;
nFrames = 0;
tic
while ~feof(x)
    data = textscan(x,sForm,nSamples,'Delimiter','\t');
    
    for ii=1:length(ste)
        sf = strfind(data{1}{1},ste{ii});
        if ~isempty(sf)
            
            switch ste{ii}
                case 'bitDepth'
                    qm = strfind(data{1}{1},'"');
                    recInfo.bitDepth = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'bitDepth')))=[];
                    break
                case 'dwellTime'
                    qm = strfind(data{1}{1},'"'); 
                    recInfo.dwellTime = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'dwellTime')))=[];
                    break
                case 'laserPower'
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.pockels = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'laserPower')))=[];
                    break
                case 'laserWavelength'
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.laserWavelength = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'laserWavelength')))=[];
                    break
                case 'micronsPerPixel'
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    mpp(1) = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    mpp(2) = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    recInfo.micronsPerPixel = mpp;
                    ste(~cellfun(@isempty,strfind(ste,'micronsPerPixel')))=[];
                    break
                case 'objectiveLens'
                    qm = strfind(data{1}{1},'"');
                    recInfo.objective = data{1}{1}(qm(3)+1:qm(4)-1);
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.objectiveMag = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.objectiveNA = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'objectiveLens')))=[];
                    break
                case 'opticalZoom'
                    qm = strfind(data{1}{1},'"');
                    recInfo.opticalZoom = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'opticalZoom')))=[];
                    break
                case 'pmtGain'
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.PMTgain_ch01 = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.PMTgain_ch02 = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'pmtGain')))=[];
                    break
                case 'ZAxis'
                    data = textscan(x,sForm,nSamples,'Delimiter','\t');
                    qm = strfind(data{1}{1},'"');
                    recInfo.Zdepth = str2double(data{1}{1}(qm(3)+1:qm(4)-1));
                    ste(~cellfun(@isempty,strfind(ste,'ZAxis')))=[];
                    break
                case 'Sequence type'
                    qm = strfind(data{1}{1},'"');
                    recInfo.timeOfRec = data{1}{1}(qm(5)+1:qm(6)-1);
                    ste(~cellfun(@isempty,strfind(ste,'Sequence type')))=[];
                    break
                case '.ome.tif'
                    nFrames = nFrames+1;
                    
            end
        end
    end
end
recInfo.nFrames = nFrames;
fclose('all')
toc