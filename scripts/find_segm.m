function segments = find_segm(x_all,y_all,I,show)
% FIND_SEGM(I)
%
% Dividing image into smaller segments using watershed ridge lines for 
% background skeleton.
% Input:
% x_all - x-axis labels [mx1]
% y_all - y-axis labels [1xn]
% I - image [mxn]
% disk_size - size of morphological structuring element
% show - if plot results
% Output:
% segments - structure with informations about segments:
% - .data - segment image
% - .x - segment positions on x-scale
% - .y - segment positions on y-scale
%
% Author: Michal.Marczyk@polsl.pl% x - x-axis labels [mx1]
% y - y-axis labels [1xn]
warning('off','images:initSize:adjustingMag');

if nargin < 2; show = false; end
if nargin < 1; error('No image to segment.'); end

%find watershed ridge lines
disk_size = 1;
bgm = water_lines(I,1,disk_size);

%select regions
segments = select_reg(x_all,y_all,I,bgm);

%crop segments
for a=1:length(segments)    
    segments(a) = auto_crop(segments(a));
end

if show
    se2 = strel('disk', 1);
    figure; 
    imshow(I,[]); hold on;
    tmp = uint8(1-bgm); tmp(tmp==1) = 255;
    himage2 = imshow(255-imdilate(255-tmp,se2));
    tmp2 = ~imdilate(bgm,se2) * 0.1; tmp2(tmp2==0) = 1;
    alpha(himage2, tmp2)
    title([num2str(length(segments)) ' regions found using Watershed ridge lines. Disk size: ' num2str(disk_size)])
    for a=1:length(segments)
        center(1) = mean(segments(a).x);
        center(2) = mean(segments(a).y);
        text(center(2)-min(y_all)+1,center(1)-min(x_all)+1,num2str(a),'HorizontalAlignment','center')
    end
end

function bgm = water_lines(I,thr_lev,disk_size)

I = 255 - I;    %invert image
se = strel('disk', disk_size);

% Clear image
Ie = imerode(I, se);%figure; imshow(Ie)
Iobr = imreconstruct(Ie, I);%figure; imshow(Iobr)
Iobrd = imdilate(Iobr, se);%figure; imshow(Iobrd)
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));%figure; imshow(Iobrcbr)
Iobrcbr = imcomplement(Iobrcbr);%figure; imshow(Iobrcbr)

% Compute Background Markers
if thr_lev==1
    bw = im2bw(Iobrcbr, graythresh(Iobrcbr));%figure; imshow(bw)
else
    bw = im2bw(Iobrcbr, thr_lev/255);
end

%Thin" the background by computing the "skeleton by influence zones".
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;  %find watershed ridge lines

function [segments,seg_size] = select_reg(x_all,y_all,I,bgm)

regions = bwlabel(1-bgm);
n_reg = max(max(regions));
segments = struct;
seg_size = zeros(1,n_reg);
for a=1:n_reg
    %define current region mask
    segments(a).mask = regions==a;
    seg_size(a) = sum(segments(a).mask(:));
    [row,col] = find(segments(a).mask);
    
    %extract image for region
    data = I;
    data(~segments(a).mask) = 255;
    
    %crop segment
    topRow = min(row); bottomRow = max(row);
    leftCol = min(col); rightCol = max(col);
    data = data(topRow:bottomRow, leftCol:rightCol);
    x = x_all(topRow:bottomRow);
    y = y_all(leftCol:rightCol);

    segments(a).x = x;
    segments(a).y = y;
    segments(a).data = data;
end

function seg_crop = auto_crop(seg)

data = 255-seg.data;
[rows,cols] = find(data);

% Get the cropping parameters
topRow = min(rows(:));
bottomRow = max(rows(:));
leftColumn = min(cols(:));
rightColumn = max(cols(:));

% Extract a cropped image from the original.
data_crop = data(topRow:bottomRow, leftColumn:rightColumn);

seg_crop.mask = seg.mask(topRow:bottomRow, leftColumn:rightColumn);
seg_crop.x = seg.x(topRow:bottomRow);
seg_crop.y = seg.y(leftColumn:rightColumn);
seg_crop.data = 255 - data_crop;