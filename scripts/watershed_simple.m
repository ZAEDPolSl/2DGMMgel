function out = watershed_simple(I,thr)

% Convert gray image to binary with given threshold
bw = im2bw(I/255,(255-thr)/255);

% Compute the distance transform of the binary image
D = bwdist(~bw);

% Complement the distance transform,and force pixels that don't belong to the objects to be at -Inf
D = -D; D(~bw) = -Inf;

% Compute the watershed transform and set background pixels to 0
out = watershed(D);
out(out==1) = 0;