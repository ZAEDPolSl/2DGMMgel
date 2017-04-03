function gmm_new = new_scale(x,y,data,gmm)

data = double(data);
[x_vec,y_vec,stats] = par_vec(x,y);
coord_vec = [x_vec',y_vec'];

gmm_new = gmm;
parfor a=1:gmm.KS
    tmp = gmm.alpha(a)*norm_pdf2(coord_vec,gmm.center(a),gmm.covar(:,:,a));
    tmp2 = reshape(tmp,stats.len_row,stats.len_col) .* data;
    scale(a) = sum(tmp2(:));
end
disp('ok')
figure; plot(scale,gmm.alpha,'.')

gmm_new.alpha = scale;

plot_gmm(x,y,data,gmm_new,default_gmm_opts());

function y = norm_pdf2(x,center,covar)
den = (6.283185307179585 * sqrt(norm(covar)));
x_centr = bsxfun(@minus,x,center);
% x_centr = (x-repmat(center,length(x),1));
y = exp(-0.5*sum((x_centr/covar).*x_centr,2))/den;
y = y';