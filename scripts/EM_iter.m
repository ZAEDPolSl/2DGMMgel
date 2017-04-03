function gmm = EM_iter(x,y,init,opts)
% EM_ITER(X,Y,INIT)
% EM algorithm for 2D Gaussian mixture model.
% Input:
% x - pixel coordinates [n,2]
% y - intensity values [n,1]
% init - vector of initial parameters for Gaussian components
% opts - structure with algorithm parameters
% Output:
% gmm - structure with estimated model parameters and statistics
% 
% Author: Michal.Marczyk@polsl.pl

if nargin<3; error('Not enough input arguments.'); end 
if nargin<4; opts = default_gmm_opts(); end
if ~opts.show
    del = y==0;
    x(del,:) = [];
    y(del) = [];
end

y = y(:)';
N = length(y);        % data length
change = Inf;
L_new = 0;
count = 1;
KS = init.KS;
sumy = sum(y);
alpha = init.alpha;
center = init.center;
covar = init.covar;
center_thr = [min(x,[],1) max(x,[],1)];   %if component center at the border
opts.SW = eye(2)*opts.SW;

if opts.show
    fprintf(1,'EM iter. no. 1. Change:infinity'); 
    f_handle = figure;
end

% MAIN LOOP
while change > opts.eps_change && count < opts.count
    
    %removing small and border components 
    ind = alpha < opts.thr_alpha | ...
        (center_thr(1) == round(center(:,1)))' | (center_thr(2) == round(center(:,2)))' | ...
        (center_thr(3) == round(center(:,1)))' | (center_thr(4) == round(center(:,2)))';
    alpha(ind) = []; center(ind,:) = []; covar(:,:,ind) = [];
    KS = length(alpha);
    L_old = L_new;

    %calculate density function
    f = zeros(KS,N);
    for a=1:KS
        f(a,:) = norm_pdf(x,center(a,:),covar(:,:,a));
    end
    
    if opts.show
        plot_gmm_short(x,y,f,alpha,center,covar,KS,f_handle,opts.cov_type);
        pause();
    end
    px = alpha * f;
    px(isnan(px) | px==0) = 5e-324;

    for a=1:KS
        pk = (alpha(a)*f(a,:).*y)./px;    %calculate posterior probability
        denom = sum(pk);
        
        alpha(a) = denom/sumy;   
        center(a,:) = (pk*x)/denom;
        
        %sphere covariance
        covarnum = (pk*(bsxfun(@minus,x,center(a,:)).^2));
        sig2_tmp = mean(covarnum/denom);
        covar(:,:,a) = opts.SW + eye(2)*sig2_tmp;
        
        %check if singularity appeared
        if det(covar(:,:,a)) <= 0.1
            covar(:,:,a) = covar(:,:,a) .* eye(2);  
            alpha(a) = 0;
        end
    end
    L_new = sum(log(px));
    change = 100*abs((L_new - L_old)/L_old);
        
    if opts.show
        tmp = ceil(log10(count)) + 18;
        for a = 1:tmp; fprintf(1,'\b'); end
        fprintf(1,'%d. Change: %8.2g',count,change);
    end
    count = count+1;
end

% RETURN RESULTS
gmm.alpha = alpha;
gmm.center = center;
gmm.covar = covar;
gmm.KS = KS;
gmm.logL = L_new;
gmm = BIC_filt(gmm,f,sumy,x);
gmm.IC = abs(2*L_new) + (7*KS-1)*log(sumy);
gmm.IC_name = 'BIC';

function y = norm_pdf(x,center,covar)
den = (6.283185307179585 * sqrt(det(covar)));
x_centr = bsxfun(@minus,x,center);
num = -0.5*sum((x_centr/covar).*x_centr,2);
% y = exp(max(num,-745))/den;
y = exp(num)/den;
y = y';

function gmm_proc = BIC_filt(gmm,f,sumy,x)
%find components to remove or merge using BIC criterion

alpha = gmm.alpha; center = gmm.center; covar = gmm.covar; KS = gmm.KS;
L_new = gmm.logL;
run = 1;
gmm_proc = gmm;

while run && KS>1
    
    BIC = abs(2*L_new) + (7*KS-1)*log(sumy);   
    BIC_tmp = zeros(1,KS); L_tmp = BIC_tmp;
    for a=1:KS  %delete each component and calculate BIC
        ind_tmp = true(1,KS); ind_tmp(a) = false;
        alpha_tmp = alpha(ind_tmp)/sum(alpha(ind_tmp));
        px_tmp = alpha(ind_tmp)*f(ind_tmp,:);
        px_tmp(isnan(px_tmp) | px_tmp==0) = 5e-324;
        L_tmp(a) = sum(log(px_tmp));
        BIC_tmp(a) = abs(2*L_tmp(a)) + (7*(KS-1)-1)*log(sumy);   
    end
    
    %merge the closest components and calculate BIC
    ind_merge = kmeans(center,KS-1,'emptyaction','singleton','replicates',20);
    for a=1:KS-1
        ind_tmp = ind_merge==a;
        if sum(ind_tmp)>1
            ind_m = find(ind_tmp);
            alpha_tmp = alpha(ind_m(1)) + alpha(ind_m(2));
            center_tmp = (center(ind_m(1),:)*alpha(ind_m(1)) + center(ind_m(2),:)*alpha(ind_m(2)))/alpha_tmp;
            covar_tmp = (covar(:,:,ind_m(1))*alpha(ind_m(1)) + covar(:,:,ind_m(2))*alpha(ind_m(2)))/alpha_tmp;
            f_tmp = norm_pdf(x,center_tmp,covar_tmp);
            alpha_tmp2 = [alpha(~ind_tmp),alpha_tmp]/sum([alpha(~ind_tmp),alpha_tmp]);
            px_tmp = alpha_tmp2*[f(~ind_tmp,:);f_tmp];
            px_tmp(isnan(px_tmp) | px_tmp==0) = 5e-324;
            L_tmp(end+1) = sum(log(px_tmp));
            BIC_tmp(end+1) = abs(2*L_tmp(end)) + (7*(KS-1)-1)*log(sumy); 
        end
    end
    
    if sum(BIC_tmp<BIC)
        if BIC_tmp(end) == min(BIC_tmp)
            alpha(ind_m(1)) = alpha_tmp;
            center(ind_m(1),:) = center_tmp;
            covar(:,:,ind_m(1)) = covar_tmp;
            alpha(ind_m(2)) = []; center(ind_m(2),:) = []; covar(:,:,ind_m(2)) = [];
            alpha  = alpha/sum(alpha);
            f(ind_m(2),:) = []; KS = KS-1;
            L_new = L_tmp(end);
        else
            [~,ind] = min(BIC_tmp); ind = ind(1);
            alpha(ind) = []; center(ind,:) = []; covar(:,:,ind) = [];
            alpha  = alpha/sum(alpha);
            f(ind,:) = []; KS = KS-1;
            L_new = L_tmp(ind);
        end
    else
        gmm_proc.alpha = alpha;
        gmm_proc.center = center;
        gmm_proc.covar = covar;
        gmm_proc.KS = KS;
        gmm_proc.logL = L_new;
        run  = 0;
    end
end

