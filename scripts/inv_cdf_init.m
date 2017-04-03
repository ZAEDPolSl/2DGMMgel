function [alpha,mu,sigma] = inv_cdf_init(x,y,K)
%Initialization of GMM parameters by inverse CDF method

mu = zeros(1,K);
sigma = zeros(1,K);
N = length(x);
cump = zeros(1,N);
for a=1:N
    cump(a)=sum(y(1:a))/sum(y); 
end
cump(1) = 0;

for a=1:K
    see = (a-0.5)/K;
    ind = find(cump<=see, 1, 'last' );
    if ind < N
        mu(a) = x(ind)+(see-cump(ind))*(x(ind+1)-x(ind));
    else
        mu(a) = x(ind);
    end
end
poz=[x(1) mu x(end)];
for a=1:K
    sigma(a) = poz(a+1) - poz(a);
end
alpha = ones(1,K)/K;    %uniform distribution  