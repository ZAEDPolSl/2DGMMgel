function plot_gmm_short(coord,ints,f_pdf,alfa,center,covar,KS,f_handle,cov_type)

% Plot 2D gel image versus GMM decomposition. 
x = unique(coord(:,1));
y = unique(coord(:,2));
m = length(x); n = length(y);
data = reshape(255-ints,n,m)';
ploty = zeros(m,n);
for a=1:KS
    tmp = alfa(a)*f_pdf(a,:);
    ploty = ploty + reshape(tmp,n,m)';
end

figure(f_handle); 
imshow(data,[],'Border','tight','InitialMagnification','fit')
hold on;
for a=1:KS
    plot_2DGauss(center(a,:)-[min(x)-1,min(y)-1],rot90(covar(:,:,a),2),cov_type);
end
% title([num2str(KS) ' components model'])


drawnow;