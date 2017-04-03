function data_filt = proc_seg(data,opts)

data = double(data);

%filter image
switch opts.filt_type
    case 'none'
        data_filt = data;
    case '2DMF'
        data_filt = filt_match2D(data,opts);
    case 'MEDF'
        data_filt = medfilt2(data,[opts.mask,opts.mask]);
    case 'MMWF'
        [~,~,~,~,data_filt] = MMWF_2D_website(data,opts.mask); 
    case 'WAVE'
        data_filt = filt_wave(data,opts);
    otherwise
        error('Wrong filtering method.')
end

switch opts.back_type
    case 'none'
 
    case 'IPF'
        data_filt = it_polyfit(data_filt,1:size(data,1),1:size(data,2),opts.poly_deg);
    case 'RB'
        se = strel('ball',opts.ball_size,1); 
        data_filt = imtophat(data_filt,se);
    case 'OTSU'
        data2 = uint8(data_filt(:));
        level = graythresh(data2(:));
        level2 = graythresh(data2(data2<level*max(data2)));
        thr = level2*max(data2);
        data_filt = data_filt - double(thr); 
        data_filt = max(data_filt,0);        
    otherwise
        error('Wrong background correction method.')
end

%IPF - iterative polynomial fitting
function [data_filt,background] = it_polyfit(data,x,y,poly_deg)
%Iterative fitting of r'th order polynomial to each row and column of the
%2D image for background and streaks removal
x = x(:); y = y(:)';
data = double(data);
iter_nb = 500;      %number of iterations
thr_perc = 0.001;
[m,n] = size(data);
background1 = zeros(m,n);
background2 = background1;
parfor a=1:m
    warning off;
    data_tmp = data(a,:);
    f1 = zeros(1,n);
    b = 1;
    thr_val = thr_perc * max(data_tmp);
    while b<iter_nb && abs(mean(data_tmp-f1)) > thr_val
        p = polyfit(y,data_tmp,poly_deg);
        f1 = polyval(p,y);
        ind = data_tmp > f1;
        data_tmp(ind) = f1(ind);
        b = b + 1;
    end
    background1(a,:) = max(f1,0);
end
parfor a=1:n
    warning off;
    f2 = zeros(m,1);
    data_tmp = data(:,a);
    b = 1;
    thr_val = 0;%0.001 * max(data_tmp);
    while b<iter_nb && abs(mean(data_tmp-f2)) > thr_val
        p = polyfit(x,data_tmp,poly_deg);
        f2 = polyval(p,x);
        ind = data_tmp > f2;
        data_tmp(ind) = f2(ind);
        b = b + 1;
     end
    background2(:,a) = max(f2,0);
end
background = (background1 + background2)/2;
data_filt = data - background; data_filt = max(data_filt,0);

%2DMF - 2D matched filtering
function data_filt = filt_match2D(data,opts)
opts.min_spread = 1;
%step1 - dividie image by intensity scale into sub-images
dtmp = data(:); dtmp(dtmp==0) = []; dtmp(1) = 0;
data2 = cell(1,opts.win_num);
thresh= [min(dtmp),multithresh(dtmp,opts.win_num-1),max(dtmp)];
for a=1:opts.win_num
    tmp_range(1) = thresh(a) - opts.win_lap*(thresh(a+1)-thresh(a))/100;
    tmp_range(2) = thresh(a+1) + opts.win_lap*(thresh(a+1)-thresh(a))/100;
    data_tmp = data;
    data_tmp(data_tmp < tmp_range(1)) = 0;
    data2{a} = match_filt(data_tmp,opts);
end
data_filt = zeros(size(data)); for a=1:opts.win_num; data_filt = data_filt+data2{a};end

function data2 = match_filt(data,opts)
max_size = 30;
PSNR = zeros(1,max_size);
data_tmp = zeros(size(data,1),size(data,2),max_size);
for b=opts.min_spread:max_size
    gfilt = fspecial('Gaussian',2*ceil(2*b)+1,b);
    gfilt = gfilt - mean(mean(gfilt));  %setting mean to 0 at z axis
    data_tmp(:,:,b) = imfilter(data,gfilt,'replicate','conv');
    PSNR(b) = measerr(data,data_tmp(:,:,b));
end
data_tmp = max(0,data_tmp);
data2 = zeros(size(data));
[~,ind] = max(PSNR); 
if ind > opts.min_spread
    data2 = data_tmp(:,:,ind);
end

%MMWF - Median modified Wiener Filter
function [x_mean,x_median,x_min,x_wiener,x_mmwf,x_mmwf_star]=MMWF_2D_website(x,n)
warning off;
% Code designed by Carlo Vittorio Cannistraci 2014
% INPUT: x is the signal as a 2D matrix, n is the window size (only odd numbers accepted)
% OUTPUT: x_mmwf is the MMWF result; x_mmwf_star is the MMWF* result
tt=numel(x);
xf_mean=zeros(tt,1);
xf_median=xf_mean;
xf_min=xf_mean;
xf_var=xf_mean;
xf_var_median1=xf_mean;
n=(n-1)/2;
[s1,s2] = ind2sub(size(x),1:tt);
parfor i=1:tt
    xt=x(max(1,s1(i)-n):min(s1(i)+n,size(x,1)),max(1,s2(i)-n):min(s2(i)+n,size(x,2)));
    xf_mean(i)=mean(xt(:));
    xf_median(i)=median(xt(:));
    xf_min(i)=min(xt(:));
    xf_var(i)=var(xt(:));
    xf_var_median1(i)=mean((xt(:)-xf_median(i)).^2);
end
x_mean=zeros(size(x));
x_mean(1:tt)=xf_mean;

x_median=zeros(size(x));
x_median(1:tt)=xf_median;

x_min=zeros(size(x));
x_min(1:tt)=xf_min;

x_wiener=zeros(size(x));
vv=max(0,xf_var-mean(xf_var))./xf_var; vv(isnan(vv))=0;
x_wiener(1:tt)=xf_mean+vv.*(x(:)-xf_mean);

x_mmwf=zeros(size(x));
x_mmwf(1:tt)=xf_median+vv.*(x(:)-xf_median);

%%% New part
x_mmwf_star=zeros(size(x));
x_mmwf_star(1:tt)=xf_median+vv2(xf_var_median1).*(x(:)-xf_median);

%%% Support function
function vv2=vv2(xf_var)
vv2=max(0,xf_var-median(xf_var))./xf_var;
vv2(isnan(vv2))=0;

%WAVE - 2D wavelet filter
function [data_filt,thr] = filt_wave(data,opts)
opts.thr_type = 'penalized';
opts.keepapp = 1;
%perform 2D wavelet decomposition
[C,S] = wavedec2(data,opts.level,opts.wname);
switch opts.thr_type
    case 'penalized'    %Birgé-Massart
        % Estimate the noise standard deviation from the
        % detail coefficients at given level .
        det1 = detcoef2('compact',C,S,opts.level);
        sigma = median(abs(det1))/0.6745;
        thr = wbmpen(C,S,sigma,opts.alpha);
    case 'default'  %Donoho and Johnstone      
        thr = ddencmp('den','wv',data);
end
data_filt = wdencmp('gbl',C,S,opts.wname,opts.level,thr,'s',opts.keepapp);
data_filt = max(0,data_filt);