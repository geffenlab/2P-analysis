% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter. 

i = 0;

% i = i+1;
% db(i).mouse_name    = 'K071';
% db(i).date          = '20170728K071_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 12; 

i = i+1;
db(i).mouse_name    = 'K070';
db(i).date          = '20170730K070_tifStacks';
db(i).expts         = [2];
db(i).diameter      = 12; 

%% good diameters
% 8 is good for 1X zoom
%  12  for 1.25X
% 14 for 1.5X
% 30 for 4X

% 
% i = i+1;
% db(i).mouse_name    = 'K058';
% db(i).date          = '20170523K058_tifStacks';
% db(i).expts         = [2];
% db(i).diameter      = 30; 

% i = i+1;
% db(i).mouse_name    = 'K056';
% db(i).date          = '20170421K056_tifStacks';
% db(i).expts         = [2];
% db(i).diameter      = 12; % 8 is good for 1X recording, 14 for 1.5X

% 
% i = i+1;
% db(i).mouse_name    = 'K057';
% db(i).date          = '20170417K057_tifStacks';
% db(i).expts         = [2];
% db(i).diameter      = 12; % 8 is good for 1X recording, 14 for 1.5X

% db(i).RootDir       = 'C:\data\K048\20170313K048_tifStacks\03132017-FMtones\'

% for ii=1:length(db(i).expts)
%     folderNames{ii} = ['C:\data\' db(i).mouse_name '\' db(i).date '\' num2str(db(i).expts(ii)) '\'];
% end

% n = getNumbersOfFrames(folderNames);

% save(['C:\data\' db(i).mouse_name '\' db(i).date '\framesPerRec.mat'],'n')


%% 
% i = i+1;
% db(i).mouse_name    = 'K048';
% db(i).date          = '20170123K048_tifStacks';
% db(i).expts         = [3];
% db(i).diameter      = 14;
% 
% 
% i = i+1;
% db(i).mouse_name    = 'K019';
% db(i).date          = '20161027K019_tifStacks'; 
% db(i).expts         = [3]; % Name of folder with tif-stack in
% db(i).nplanes       = 1; 
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% % db(i).expred        = [3];
% 
% 
% i = i+1;
% db(i).mouse_name    = 'K019';
% db(i).date          = '20161027K019_tifStacks'; 
% db(i).expts         = [4]; % Name of folder with tif-stack in
% db(i).nplanes       = 1; 
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% % db(i).expred        = [3];
% 
% % i = i+1;
% db(i).mouse_name    = 'M150329_MP009';
% db(i).date          = '2015-04-29';
% db(i).expts         = [4 5 6];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1; 
% db(i).expred        = [3];
% db(i).nchannels_red = 2;
% db(i).comments      = 'multi p file: block 4,5,6';
% 
% i = i+1;
% db(i).mouse_name    = 'M150331_MP011';
% db(i).date          = '2015-04-29';
% db(i).expts         = [3 4 6];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1; 
% db(i).expred        = [2];
% db(i).nchannels_red = 2;
% db(i).comments      = 'multi p file: block 4,3,6';
% 
% i = i+1;
% db(i).mouse_name    = 'M150422_MP012';
% db(i).date          = '2015-04-28';
% db(i).expts         = [3 4 5];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1; 
% db(i).expred        = [2];
% db(i).nchannels_red = 2;
% db(i).comments      = 'multi p file: block 3,4,5';
% 
% i = i+1;
% db(i).mouse_name    = 'M150422_MP012';
% db(i).date          = '2015-05-04';
% db(i).expts         = [2 8];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1; 
% db(i).comments      = 'single p file: block 2';
% 
% i = i+1;
% db(i).mouse_name    = 'M150422_MP012';
% db(i).date          = '2015-05-20';
% db(i).expts         = [4];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1; 
% db(i).expred        = [2];
% db(i).nchannels_red = 2;
% db(i).comments      = 'single p file: block 4';
% 
% i = i+1;
% db(i).mouse_name    = 'M150422_MP015';
% db(i).date          = '2015-05-09';
% db(i).expts         = [2 3];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1;
% db(i).expred        = [1];
% db(i).nchannels_red = 2;
% db(i).comments      = 'single p file, block 3';
% 
% % i = i+1;
% % db(i).mouse_name    = 'M150808_MP016';
% % db(i).date          = '2015-08-24';
% % db(i).expts         = [4];
% % db(i).nchannels     = 1;
% % db(i).gchannel      = 1; 
% % db(i).nplanes       = 1; 
% % db(i).expred        = [3];
% % db(i).nchannels_red = 2;
% % db(i).comments      = 'single p file, block 4, zoom 20';
% 
% i = i+1;
% db(i).mouse_name    = 'M150423_MP014';
% db(i).date          = '2015-06-16';
% db(i).expts         = [8];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1; 
% db(i).expred        = [4];
% db(i).nchannels_red = 2;
% db(i).comments      = 'single p file, block 8 + .3Hz stimulus';
% 
% i = i+1;
% db(i).mouse_name    = 'M150422_MP015';
% db(i).date          = '2015-05-01';
% db(i).expts         = [5 6];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1;
% db(i).expred        = [];
% db(i).nchannels_red = 2;
% db(i).comments      = 'only 1Hz vs 10Hz';
% 
% i = i+1;
% db(i).mouse_name    = 'M150331_MP011';
% db(i).date          = '2015-05-02';
% db(i).expts         = [4 7 8 10 11 12];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1;
% db(i).expred        = [3];
% db(i).nchannels_red = 2;
% db(i).comments      = 'first three running, last three not running';
% 
% i = i+1;
% db(i).mouse_name    = 'M150422_MP015';
% db(i).date          = '2015-04-28';
% db(i).expts         = [3 4 5];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1;
% db(i).expred        = [];
% db(i).nchannels_red = 2;
% db(i).comments      = 'might not be enough signal in this one!';
% 
% i = i+1;
% db(i).mouse_name    = 'M150808_MP016';
% db(i).date          = '2015-08-24';
% db(i).expts         = [4];
% db(i).nchannels     = 1;
% db(i).gchannel      = 1; 
% db(i).nplanes       = 1;
% db(i).expred        = 3;
% db(i).nchannels_red = 2;
% db(i).comments      = 'zoom 20x recording';
