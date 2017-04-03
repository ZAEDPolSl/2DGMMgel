function [x,ind] = Dixon_test(x,sig_lev)
% DIXON_TEST(X,SIG_LEV)
% Inputs:
% x - vector of data [1xn]
% sig_lev - significance level
% Output:
% x - vector of cleaned data [1xn or 1xn-1]
% ind - index of outlier, if any
%
% Author: Michal.Marczyk@polsl.pl

if nargin<2; sig_lev = 0.05; end
if nargin<1; error('Not enough input arguments'); end

[x,ind2] = sort(x);

if sig_lev == 0.1
    Qcrit = [1,1,0.941,0.765,0.642,0.560,0.507,0.468,0.437,0.412];
elseif sig_lev == 0.05
    Qcrit = [1,1,0.970,0.829,0.710,0.625,0.568,0.526,0.493,0.466];
elseif sig_lev == 0.01
    Qcrit = [1,1,0.994,0.926,0.821,0.740,0.680,0.634,0.598,0.568];
end
	
ind = false(1,length(x));

r10(1) = (x(2)-x(1))/(x(end)-x(1));
r10(2) = (x(end)-x(end-1))/(x(end)-x(1));

% removing only one border outlier
if  r10(1) > r10(2)
    if r10(1) > Qcrit(length(x))
        ind(ind2(1)) = true;
    end
else
    if  r10(2) > Qcrit(length(x))
        ind(ind2(end)) = true;
    end
end
x(ind) = [];
