function f = fast_norm_pdf(x,center,covar)

N = size(x,1);
KS = size(center,1);

f = zeros(N,KS);
for a=1:KS
    den = 6.283185307179585 * sqrt(norm(covar(:,:,a)));
    x_centr = bsxfun(@minus,x,center(a,:));
    f(:,a) =  exp(-0.5*sum((x_centr/covar(:,:,a)).*x_centr,2))/den;
end
f = f';
