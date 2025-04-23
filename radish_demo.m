%% Example run to generate B1 and B0 maps from WASABI data

addpath(genpath('./')); % replace with path to RADISH

%% Load and normalise your WASABI data 
% load_nii requires the following package by Jimmy Shen: https://au.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image
wasabi = load_nii('wasabi.nii.gz'); % replace with link to your WASABI data
mask = load_nii('mask.nii.gz');
% or: load('wasabi.mat');

% If 3D data x,y,z,t
wnorm = double(wasabi.img(:,:,:,2:end)) ./ double(wasabi.img(:,:,:,1));
% else if 2D data x,y,t
% wnorm = double(wasabi.img(:,:,2:end)) ./ double(wasabi.img(:,:,1));

%% Initialise supporting struct
w = setup_radish;

% modify the below parameters to suit your data
w.tp = 0.00512; % s
w.w0 = 298; % MHz
w.max = 1.5; % maximum offset value [ppm]
w.noffsets = 49; % number of offsets

%% Run RADISH 
[b1, b0, scores, errored] = radish(wnorm, mask, w);

rb1 = b1 / 3.7; % change to nominal B1 value if needed

%% Quality control
qa_thresh = 10; % set the accepted threshold for score
qa_failed = zeros(size(b1));
qa_failed(scores > qa_thresh) = 1;
% Try rerunning, this time ignoring B0 values that are problematic
excludeB0 = b0(qa_failed == 1);
[b1_rerun, b0_rerun, scores_rerun, errored_rerun] = radish(wnorm, mask, w, excludeB0);
