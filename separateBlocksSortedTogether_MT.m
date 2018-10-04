% separate blocks sorted together:

% MT Modification 2018/04/19 : can move the exptInfo files correctly if
% there are several blocks for a same recording
% Note : the recInfo read from xml file will only go to the first exptInfo
% file (block01)

clear
mouse = {'MT004'};
date = '20180417';

for zz = 1:length(mouse)
    fold_loc = ['E:\dataAnalysed\' mouse{zz} '\' date mouse{zz} '_tifStacks\'];
    folders = dir(fold_loc);
    folders(strcmp({folders.name},'.'))=[]; folders(strcmp({folders.name},'..'))=[];
    folders(~contains({folders.name},'_'))=[]; % only look at folders with underscores!
    
    for uu = 1:length(folders)
        % define folders:
        us = strfind(folders(uu).name,'_');
        f{1} = folders(uu).name(1:us(1)-1);
        for ii=2:length(us)
            f{ii} = folders(uu).name(us(ii-1)+1:us(ii)-1); %#ok<SAGROW>
        end
        f{end+1} = folders(uu).name(us(end)+1:end); %#ok<SAGROW>
        
        %% Move the xml files ...
        for ii=1:length(f)
            fLoc = ['C:\data\' mouse{zz} '\' '\' date mouse{zz} '_tifStacks\' f{ii} '\'];
            if ~isdir([fold_loc f{ii} '\']); mkdir([fold_loc f{ii} '\']); end
            files = dir(fLoc);
            files(~cellfun(@isempty,strfind({files.name},'.tif')))=[];
            files(~cellfun(@isempty,strfind({files.name},'.db')))=[];
            files([files.isdir])=[];
            for mm = 1:length(files)
                movefile([fLoc files(mm).name],[fold_loc f{ii} '\']);
            end
        end
        
       % get xml data
        for ii=1:length(f)
            xml = dir([fold_loc f{ii} '\*xml']);
            recInfo = readXmlFile_v2_20170730([fold_loc f{ii} '\' xml.name]);
            save([fold_loc f{ii} '\recInfo.mat'],'recInfo')
            ri(ii) = recInfo; %#ok<SAGROW>
        end
        
        tr = zeros(1,length(ri));
        for ii = 1:length(ri)
            colon = strfind(ri(ii).date,':');
            space = strfind(ri(ii).date,' ');
%             tr(ii) = str2double([ri(ii).date(space+1:colon(1)-1),ri(ii).date(colon(1)+1:colon(2)-1)]);
            tr(ii) = datenum(ri(ii).date,'yyyy-mm-dd HH:MM:SS');
        end
        %% move exptInfo into correct folders      
        
%         exptInfoFold = [fold_loc folders(uu).name '\*exptInfo.mat'];
%         ei = dir(exptInfoFold);
%         if ~isempty(ei)
%             t = zeros(1,length(ei));
%             for ii=1:length(ei)
%                 us = strfind(ei(ii).name,'_');
%                 t(ii) = str2double(ei(ii).name(us(1)+1:us(1)+4)); % times of expt files
%             end
%         else



            exptInfoFold = ['C:\data\' mouse{zz} '\' date mouse{zz} '_tifStacks\*exptInfo.mat'];
            exf = dir(exptInfoFold);
            t = zeros(1,length(exf));
            for ii=1:length(t)
                us = strfind(exf(ii).name,'_');
%                 t(ii) = str2double(exf(ii).name(us(1)+1:us(1)+4)); % times of expt files
                t(ii) = datenum(exf(ii).name(1:us(2)-1),'yymmdd_HHMMSS');
            end
            
            tm = zeros(1,length(tr));
            for ii = 1:length(tr)
                x = etime(datevec(t),datevec(tr(ii)));
                tm(ii) = knnsearch(x,0);
%                 tm(ii) = knnsearch(t',tr(ii));
                movefile(['C:\data\' mouse{zz} '\' date mouse{zz} '_tifStacks\' exf(tm(ii)).name], [fold_loc f{ii} '\'])
            end
    
            

                for jj=1:length(tr)-1
                    [sorttr indtr] = sort(tr); % find chronological order
                   exf = dir(exptInfoFold); 
                   t = zeros(1,length(exf));
                   for ii=1:length(t)
                     us = strfind(exf(ii).name,'_');
                     t(ii) = datenum(exf(ii).name(1:us(2)-1),'yymmdd_HHMMSS');
                   end
                a=find((t>tr(indtr(jj))) & (t<tr(indtr(jj+1))));
                for ii=1:length(a)
                movefile(['C:\data\' mouse{zz} '\' date mouse{zz} '_tifStacks\' exf(a(ii)).name], [fold_loc f{indtr(jj)} '\']) 
                end
                end
                for jj=length(tr)
                   exf = dir(exptInfoFold); 
                   t = zeros(1,length(exf));
                   for ii=1:length(t)
                     us = strfind(exf(ii).name,'_');
                     t(ii) = datenum(exf(ii).name(1:us(2)-1),'yymmdd_HHMMSS');
                   end
                a=find((t>tr(indtr(jj))));
                for ii=1:length(a)
                movefile(['C:\data\' mouse{zz} '\' date mouse{zz} '_tifStacks\' exf(a(ii)).name], [fold_loc f{indtr(jj)} '\']) 
                end
                end
                        
%                 if (t(ii)>tr(end))
%                 movefile(['C:\data\' mouse{zz} '\' date mouse{zz} '_tifStacks\' exf(ii).name], [fold_loc f{length(tr)} '\']) 
%                 end
%         end
      
             
        proc = dir([fold_loc folders(uu).name '\*_proc.mat']);
        load([fold_loc folders(uu).name '\' proc.name]);
        nframes = [0,dat.ops.Nframes];
        bf = dat.ops.badframes;
        badFrames = cell(1,length(nframes)-1);
        for ii = 2:length(nframes)
            badFrames{ii-1} = bf(sum(nframes(1:ii-1))+1:sum(nframes(1:ii)));
        end
            
        datAll = dat;
        [~, b]=sort(indtr);
        for qq = 1:length(f)
            ii = b(qq);
            dat = datAll; 
            dat.Fcell = dat.Fcell(ii); 
            dat.FcellNeu = dat.FcellNeu(ii); 
            dat.sp = dat.sp{ii}; 
            dat.ca = dat.ca{ii}; 
            dat.kernels = dat.kernels{ii};
            dat.badFrames = badFrames{ii};
            save([fold_loc f{ii} '\' proc.name],'dat');
        end
        
   
        
        %% Add recInfo to exptInfo
        
        for ii=1:length(f)
            eit = dir([fold_loc f{ii} '\*exptInfo.mat']);
            load([fold_loc f{ii} '\' eit(1).name])
            rci = dir([fold_loc f{ii} '\recInfo.mat']);
            load([fold_loc f{ii} '\' rci.name])
            exptInfo.recInfo = recInfo;
            delete([fold_loc f{ii} '\' rci.name])
            save([fold_loc f{ii} '\' eit(1).name],'exptInfo','-append')
        end
        
        
        
    end
    
    clearvars -except mouse date
end

disp('Finished separating blocks');






