% combine some blocks
clear

mouse = 'K056';
date = '20170425';
folders = {'2','3'};
fold_loc = ['E:\dataAnalysed\' mouse '\' date mouse '_tifStacks\'];
nf = [fold_loc sprintf('%s_%s',folders{1},folders{2})];
mkdir(nf)
for ii=1:length(folders)
    ei = dir([fold_loc folders{ii} '\*exptInfo.mat']);
    a(ii) = load([fold_loc folders{ii} '\' ei.name]);
    proc = dir([fold_loc folders{ii} '\*proc.mat']);
    p(ii) = load([fold_loc folders{ii} '\' proc.name]);
end

fs = a(1).fs;
frameOn = [a(1).frameOn(1:length(p(1).dat.Fcell{1})); a(2).frameOn + a(1).frameOn(length(p(1).dat.Fcell{1}))];
frameOff = [a(1).frameOff(1:length(p(1).dat.Fcell{1})); a(2).frameOff(1:end) + a(1).frameOff(length(p(1).dat.Fcell{1}))];
exptInfo = a(1).exptInfo(1);
eventOff = [a(1).eventOff; a(2).eventOff(2:end) + a(1).eventOff(end)];
eventOn = [a(1).eventOn; a(2).eventOn + a(1).eventOn(end)];

save([nf 'combined_' ei.name],'exptInfo','fs','frameOn','frameOff','eventOn','eventOff')

dat = p(1).dat;
dat.Fcell{1} = [p(1).dat.Fcell{1},p(2).dat.Fcell{1}];
dat.FcellNeu{1} = [p(1).dat.FcellNeu{1},p(2).dat.FcellNeu{1}];
save([nf 'combined_' proc.name],'dat')
