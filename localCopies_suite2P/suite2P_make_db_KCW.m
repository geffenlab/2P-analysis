% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter. 

% ******* good diameters  *********
% 8 is good for 1X zoom
% 12  for 1.25X
% 14 for 1.5X
% 30 for 4X

i = 0;

%% Melanie's
% % % % % % 
% i = i+1;
% db(i).mouse_name    = 'MT019';
% db(i).date          = '20180918MT019_tifStacks';
% db(i).expts         = [1 2 3 4 5 6 7 8];
% db(i).diameter      = 12;
% % % % % % % % % 
% i = i+1;
% db(i).mouse_name    = 'MT026';
% db(i).date          = '20180928MT026_tifStacks';
% db(i).expts         = [2 3 4];
% db(i).diameter      = 12;
% % % % % % % 
% i = i+1;
% db(i).mouse_name    = 'MT017';
% db(i).date          = '20180808MT017_tifStacks';
% db(i).expts         = [1];
% db(i).expred        = [1]; % which experiments are relevant for the red
% db(i).nchannels_red = 1; % which expt (folder) is the red channel
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 14;

%% Kath's
% % 
% % % 
i = i+1;
db(i).mouse_name    = 'K133';
db(i).date          = '20180917K133_tifStacks';
db(i).expts         = [1 2];
db(i).diameter      = 12;
% % 
% % 
% i = i+1;
% db(i).mouse_name    = 'K134';
% db(i).date          = '20181003K134_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 14;
% % 
% i = i+1;
% db(i).mouse_name    = 'K134';
% db(i).date          = '20180919K134_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 14;
% 
% i = i+1;
% db(i).mouse_name    = 'K118';
% db(i).date          = '20180617K118_tifStacks';
% db(i).expts         = [1 2 4];
% db(i).expred        = [1 2]; % which experiments are relevant for the red
% db(i).nchannels_red = 4; % which expt (folder) is the red channel
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel


%% Aaron's
% % % 

% i = i+1;
% db(i).mouse_name    = 'AW017';
% db(i).date          = '20180530AW017_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 14;

