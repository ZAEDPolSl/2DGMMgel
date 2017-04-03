% Modeling 2D gel images with mixture of 2D Gaussian functions by modified
% EM algorithm. Algorihtm is composed of dividing image into smaller
% segments, spot detection within segments, fitting Gaussian mixture model
% with circular shape each segment and full model post-processing.
% Author: Michal.Marczyk@polsl.pl

addpath('scripts')

%input dataset
name = 'RG1';   %name of the dataset to load
load(['data_' name '_filt.mat'])
x = 1:size(data,1);
y = 1:size(data,2);

%set GMM options
opts = default_gmm_opts();  %script with all parameters
opts.parallel = true;       %if use parallel computing
opts.show = true;           %if present results on the screen
opts.seg = true;            %if perform image fragmentation
opts.proc_type = 'global';  %local or global processing
opts.back_type = 'OTSU';    %background correction method {'none','IPF','RB','OTSU'}
opts.filt_type = 'MMWF';    %filltering method {'none','MEDF','MMWF','2DMF','WAVE'}
opts.post = true;           %if post-procesing is on
opts.init_type = 'watershed';%'watershed' or 'external' spot detection for IC of model parameters
opts.water_min = 250;       %watershed thresholding parameter
%if external method is used spot centers must be provided as
%opts.init_spots with size nx2, where n is the number of spots. In the first
%column there is a location of each spot on x-axis and in the second column
%a location of each spot on y-axis.

%run 2DGMM_gel
[gmm,spot_det] = proc_gmm(uint8(data),x,y,opts);