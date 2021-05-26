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
% db(i).mouse_name    = 'MT073';
% db(i).date          = '20200316MT073_tifStacks';
% db(i).expts         = [8];
% db(i).expred        = [8]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 12;

% i = i+1;
% db(i).mouse_name    = 'MT082';
% db(i).date          = '20210107MT082_tifStacks';
% db(i).expts         = [2];
% db(i).expred        = [2]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 14;
% i = i+1;
% db(i).mouse_name    = 'MT084';
% db(i).date          = '20210107MT084_tifStacks';
% db(i).expts         = [3];
% db(i).expred        = [3]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 14;
%  i = i+1;
%  db(i).mouse_name    = 'MT089';
%  db(i).date          = '20210310MT089_tifStacks';
%  db(i).expts         = [5];
%  db(i).diameter      = 14;

% % % % % % % % % % % % 
%% Xiaomao
% i = i+1;
% db(i).mouse_name    = 'XD028';
% db(i).date          = '20200822XD028_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 12;

% i = i+1;
% db(i).mouse_name    = 'XD030';
% db(i).date          = '20200828XD030_tifStacks';
% db(i).expts         = [1 2 3];
% db(i).diameter      = 12;
% 
% i = i+1;
% db(i).mouse_name    = 'XD028';
% db(i).date          = '20200914XD028_tifStacks';
% db(i).expts         = [1 2 3];
% db(i).diameter      = 12;
% 
% i = i+1;
% db(i).mouse_name    = 'XD030';
% db(i).date          = '20200925XD030_tifStacks';
% db(i).expts         = [1 2 3];
% db(i).diameter      = 12;

% i = i+1;
% db(i).mouse_name    = 'XD034';
% db(i).date          = '20210306XD034_tifStacks';
% db(i).expts         = [1 2 3];
% db(i).diameter      = 12;
% 
i = i+1;
db(i).mouse_name    = 'XD036';
db(i).date          = '20210422XD036_tifStacks';
db(i).expts         = [1];
db(i).diameter      = 12;

%% Kath's
% % % 
% % 
% i = i+1;
% db(i).mouse_name    = 'K169';
% db(i).date          = '20201108K169_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 12;
% 
% i = i+1;
% db(i).mouse_name    = 'K170';
% db(i).date          = '20201108K170_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 12;

% 
% 
% i = i+1;
% db(i).mouse_name    = 'K172';
% db(i).date          = '20201111K172_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 12;




%% Aaron's
% % % 
% % % % 
% i = i+1;
% db(i).mouse_name    = 'AW075';
% db(i).date          = '20190617AW075_tifStacks';
% db(i).expts         = [1 2];
% db(i).diameter      = 12;

% i = i+1;
% db(i).mouse_name    = 'AW069';
% db(i).date          = '20190427AW069_tifStacks';
% db(i).expts         = [2 3];
% db(i).diameter      = 12;
% 
% i = i+1;
% db(i).mouse_name    = 'AW071';
% db(i).date          = '20190427AW071_tifStacks';
% db(i).expts         = [2];
% db(i).diameter      = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i = i+1;
% db(i).mouse_name    = 'AW055';
% db(i).date          = '20190404AW055_tifStacks';
% db(i).expts         = [3];
% db(i).expred        = [3]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 12;

%% Alex's
% % % 
% % % % 
% i = i+1;
% db(i).mouse_name    = 'AL043';
% db(i).date          = '20200814AL043_tifStacks';
% db(i).expts         = [4];
% db(i).expred        = [4]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 8;
% 
% i = i+1;
% db(i).mouse_name    = 'AL043';
% db(i).date          = '20200814AL043_tifStacks';
% db(i).expts         = [5];
% db(i).expred        = [5]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% i = i+1;
% db(i).mouse_name    = 'AL043';
% db(i).date          = '20200825AL043_tifStacks';
% db(i).expts         = [5];
% db(i).expred        = [5]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 8;
% 
% i = i+1;
% db(i).mouse_name    = 'AL043';
% db(i).date          = '20200825AL043_tifStacks';
% db(i).expts         = [6];
% db(i).expred        = [6]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 8;
% 
% i = i+1;
% db(i).mouse_name    = 'AL043';
% db(i).date          = '20200825AL043_tifStacks';
% db(i).expts         = [7];
% db(i).expred        = [7]; % which experiments are relevant for the red
% db(i).nchannels_red = 2; % number of channels
% db(i).nchannels     = 2; % this will always be 2 if there is a red channel
% db(i).diameter      = 8;

