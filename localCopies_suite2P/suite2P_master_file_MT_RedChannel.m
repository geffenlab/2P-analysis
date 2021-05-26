%% SET ALL DEFAULT OPTIONS HERE

clear all

% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter

% check out the README file for detailed instructions (and extra options)
addpath('C:\Users\2P-analysis\Documents\GitHub\2P-analysis\localCopies_suite2P') % add the path to your make_db file
addpath(genpath('C:\Users\2P-analysis\Documents\GitHub\OASIS_matlab\')); % add deconvolution path

% overwrite any of these default options in your make_db file for individual experiments
suite2P_make_db_KCW; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE

ops0.toolbox_path = 'C:\Users\2P-analysis\Documents\GitHub\Suite2P';
if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 1; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).
ops0.fig                    = 1; % turn off figure generation with 0

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = 'C:\data\'; % Suite2P assumes a folder structure, check out README file
ops0.temp_tiff              = 'C:\_tempreg\temp.tif'; % copies each remote tiff locally first, into this file
ops0.RegFileRoot            = 'C:\suite2PtempFolder_binary_files\';  % location for binary file
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = 'E:\dataAnalysed\'; % a folder structure is created inside
ops0.RegFileTiffLocation    = ['E:\dataAnalysed\']; %'E:\dataAnalysed\'; % leave empty to NOT save registered tiffs (slow)
ops0.regFileBinLocation     = [];%'E:\dataAnalysed\'; % save binary file:?  Added by Kath 20170827

% registration options
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 
ops0.nimgbegend             = 250; % frames to average at beginning and end of blocks
ops0.dobidi                 = 1; % infer and apply bidirectional phase offset

% cell detection options
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters 
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters % usually 0.5
ops0.nSVDforROI             = 1000; % how many SVD components for cell clustering % usually 1000
ops0.NavgFramesSVD          = 5000; % how many (binned) timepoints to do the SVD based on
ops0.signalExtraction       = 'raw'; % how to extract ROI and neuropil signals: 
%  'raw' (no cell overlaps), 'regression' (allows cell overlaps), 
%  'surround' (no cell overlaps, surround neuropil model)
ops0.refine                 = 1; % whether or not to refine ROIs (refinement uses unsmoothed PCs to compute masks)

% ----- neuropil options (if 'surround' option) ------------------- %
% all are in measurements of pixels
ops0.innerNeuropil  = 1; % padding around cell to exclude from neuropil
ops0.outerNeuropil  = Inf; % radius of neuropil surround
% if infinity, then neuropil surround radius is a function of cell size
if isinf(ops0.outerNeuropil)
    ops0.minNeuropilPixels = 400; % minimum number of pixels in neuropil surround
    ops0.ratioNeuropil     = 5; % ratio btw neuropil radius and cell radius
    % radius of surround neuropil = ops0.ratioNeuropil * (radius of cell)
end

% spike deconvolution options
ops0.imageRate              = 30;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 2; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = 1; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
% ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
% ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)

% ----- if you have a RED channel ---------------------- ------------%
ops0.AlignToRedChannel      = 0; % compute registration offsets using red channel
ops0.REDbinary              = 1; % make a binary file of registered red frames
% if db.expred, then compute mean red image for green experiments with red
% channel available while doing registration
ops0.redMeanImg             = 0; 
% for red cell detection (identify_redcells_sourcery.m)
% redratio = red pixels inside / red pixels outside
% redcell = redratio > mean(redratio) + redthres*std(redratio)
% notred = redratio < mean(redratio) + redmax*std(redratio)
ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

ops0.redthres               = 1.5; % the higher the thres the less red cells
ops0.redmax                 = 1; % the higher the max the more NON-red cells

%% RUN THE PIPELINE HERE
db0 = db;

for iexp = 1:length(db) %[3:length(db) 1:2]
    db = db0(iexp);
    run_pipeline(db, ops0);
    
    % deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
    add_deconvolution(ops0, db0(iexp));
    
    % add red channel information (if it exists)
     if isfield(db,'expred') && ~isempty(db.expred)
        % creates mean red channel image aligned to green channel
        % use this if you didn't get red channel during registration
        % OR you have a separate experiment with red and green just for this
        red_expts = ismember(db.expts, getOr(db, 'expred', []));
        if ~ops0.redMeanImg || sum(red_expts)==0
            run_REDaddon_sourcery(db, ops0);
        end
        
        % identify red cells in mean red channel image
        % fills dat.stat.redcell, dat.stat.notred, dat.stat.redprob
        identify_redcells_sourcery(db, ops0); 
        
    end
    
end
%% STRUCTURE OF RESULTS FILE
% cell traces are in dat.Fcell
% neuropil traces are in dat.FcellNeu
% manual, GUI overwritten "iscell" labels are in dat.cl.iscell
%  
% stat(icell) contains all other information:
% iscell: automated label, based on anatomy
% neuropilCoefficient: neuropil subtraction coefficient, based on maximizing the skewness of the corrected trace (ICA)
% st: are the deconvolved spike times (in frames)
% c:  are the deconvolved amplitudes
% kernel: is the estimated kernel
%% stuff to do after processing
% moveAnalysedFiles
disp('finished analysing')