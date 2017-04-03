function init = watershed_init(data,x,y,opts)
% Initialization of parameters of 2D Gaussian components for EM algorithm
% using watershed algorithm

data2 = watershed_simple(data,opts.water_min);

%find spot properties
spots = unique(data2); spots(spots==0) = [];
n_spots = length(spots);
init.alpha = zeros(1,n_spots);
init.center = zeros(n_spots,2);
init.covar = zeros(2,2,n_spots);
alpha_den = length(data2(data2>0));
x_min = min(x); y_min = min(y);
for a=1:n_spots
    ind = data2==spots(a);
    init.alpha(a) = sum(sum(ind))/alpha_den;
    tmp_data = zeros(size(data2));
    tmp_data(ind) = 1;
    stats = regionprops(tmp_data,'Centroid');
    init.center(a,:) = [x_min + stats.Centroid(2)-1,y_min+ stats.Centroid(1)-1];
    init.covar(:,:,a) = diag([(max(sum(tmp_data,1))/6)^2,(max(sum(tmp_data,2))/6)^2]);
end
init.KS = n_spots;