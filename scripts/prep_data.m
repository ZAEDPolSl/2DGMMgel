function [coord,ints] = prep_data(im,x,y)
%Prepare gray image data for GMM_2DGel

if nargin < 2
    x = 1:size(im,1);
    y = 1:size(im,2);
end
N = length(im(:));

coord = zeros(N,2);
ints = zeros(N,1);
count = 1;
for a=1:length(x)
    for b=1:length(y)
        coord(count,1) = x(a);
        coord(count,2) = y(b);
        ints(count) = im(a,b);
        count = count + 1;
    end
end


