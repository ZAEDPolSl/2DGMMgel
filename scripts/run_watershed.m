function data_out = run_watershed(data,thr_lev,disk_size)
% Marker-Controlled Watershed Segmentation
% 1. Compute a segmentation function.  This is an image whose dark
% regions are the objects you are trying to segment.
% 2. Compute foreground markers.  These are connected blobs of pixels
% within each of the objects.
% 3. Compute background markers.  These are pixels that are not part of
% any object.
% 4. Modify the segmentation function so that it only has minima at the
% foreground and background marker locations.
% 5. Compute the watershed transform of the modified segmentation function.

data = 255 - data;    %invert image
se = strel('disk', disk_size);

% Clear image
Ie = imerode(data, se);%figure; imshow(Ie)
Iobr = imreconstruct(Ie, data);%figure; imshow(Iobr)
Iobrd = imdilate(Iobr, se);%figure; imshow(Iobrd)
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));%figure; imshow(Iobrcbr)
Iobrcbr = imcomplement(Iobrcbr);%figure; imshow(Iobrcbr)

%Compute Foreground Markers
fgm = imregionalmax(Iobrcbr);
se2 = strel(ones(1,1)); % (5,5)
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 1); % 20

% Compute Background Markers
bw = im2bw(Iobrcbr/255, thr_lev*max(Iobrcbr(:))/255);figure; imshow(bw)

% Thin" the background by computing the "skeleton by influence zones".
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;  %find watershed ridge lines

% Use the Gradient Magnitude as the Segmentation Function
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(data), hy, 'replicate');
Ix = imfilter(double(data), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

% Compute the Watershed Transform of the Segmentation Function.
gradmag2 = imimposemin(gradmag, bgm | fgm4);
L = watershed(gradmag2,8);

%%
% This visualization illustrates how the locations of the foreground and
% background markers affect the result.  In a couple of locations,
% partially occluded darker objects were merged with their brighter
% neighbor objects because the occluded objects did not have foreground
% markers.
%
% Another useful visualization technique is to display the label matrix
% as a color image.  Label matrices, such as those produced by
% |watershed| and |bwlabel|, can be converted to truecolor images for
% visualization purposes by using |label2rgb|.
% data = 255 - data;
% 
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure, imshow(Lrgb)
% title('Colored watershed label matrix (Lrgb)')

data_out = Lrgb;