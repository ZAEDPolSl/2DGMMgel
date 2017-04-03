function opts = default_gmm_opts()
% Setting parameters for modeling of 2D gel images.

%Global parameters
opts.show = false;      %if display and plot results
opts.parallel = true;   %if parallel computing used
opts.seg = true;        %if perform image fragmentation
opts.post = true;       %if post-processing components

%Image processing parameters
opts.proc_type = 'local';%local or global processing
opts.back_type = 'none';%background correction method {'none','IPF','RB','OTSU'}
opts.filt_type = 'none';%filltering method {'none','MEDF','MMWF','WAVE'}
opts.poly_deg = 4;      %IPF - degree of fitting polynomial
opts.ball_size = 20;    %RB - size of rolling ball element
opts.mask = 3;          %MEDF,MMWF - moving window size
opts.win_num = 5;       %2DMF - no. of overlapping fragments
opts.win_lap = 10;      %2DMF - overlap [%]
opts.alpha = 5;         %WAVE - penalization parameter
opts.wname = 'db3';     %WAVE - wavelet type
opts.level = 2;         %WAVE - decomposition level

%Spot detection parameters
opts.init_type = 'watershed';    %initial conditions method for EM ('watershed' or 'external')
opts.water_min = 250;   %watershed thresholding parameter

%EM for GMM parameters
opts.eps_change = 1e-4; %stop criterion threshold
opts.count = 5000;      %max. no. of iterations
opts.SW = 0.1;          %regularizing coefficient for covariance
opts.thr_alpha = 0;     %small alpha threshold

%Post-processing parameters
opts.whisker = 1.5;     %1.5 - mild outliers, 3 - extreme outliers