% separate blocks sorted together:
clear
mouse = {'MT021'};
date = '20180821';

for zz = 1:length(mouse)
    fold_loc = ['E:\dataAnalysed\' mouse{zz} '\' date mouse{zz} '_tifStacks\'];
    folders = dir(fold_loc);
    folders(strcmp({folders.name},'.'))=[]; folders(strcmp({folders.name},'..'))=[];
    
    for uu = 1:length(folders)
        
        xml = dir([fold_loc folders(uu).name '\*xml']);
        recInfo = readXmlFile_v2_20170730([fold_loc folders(uu).name '\' xml.name]);
        save([fold_loc folders(uu).name '\recInfo.mat'],'recInfo')
        ri = recInfo;
        
        
        tr = zeros(1,length(ri));
        for ii = 1:length(ri)
            colon = strfind(ri(ii).date,':');
            space = strfind(ri(ii).date,' ');
            %             tr(ii) = str2double([ri(ii).date(space+1:colon(1)-1),ri(ii).date(colon(1)+1:colon(2)-1)]);
            tr(ii) = datenum(ri(ii).date,'yyyy-mm-dd HH:MM:SS');
        end
        %% move exptInfo into correct folders
        
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
            movefile(['C:\data\' mouse{zz} '\' date mouse{zz} '_tifStacks\' exf(tm(ii)).name], [fold_loc folders(uu).name '\'])
        end
        
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
        for qq=1
            ii = tm(qq);
            dat = datAll;
            dat.badFrames = badFrames;
            save([fold_loc folders(uu).name '\' proc.name],'dat');
        end
        
        
        
        %% Add recInfo to exptInfo
        
        eit = dir([fold_loc folders(uu).name '\*exptInfo.mat']);
        load([fold_loc folders(uu).name '\' eit.name])
        rci = dir([fold_loc folders(uu).name '\recInfo.mat']);
        load([fold_loc folders(uu).name '\' rci.name])
        exptInfo.recInfo = recInfo;
        delete([fold_loc folders(uu).name '\' rci.name])
        save([fold_loc folders(uu).name '\' eit.name],'exptInfo','-append')
        
        
    end
    
    clearvars -except mouse date
end








