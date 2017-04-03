function seg_crop = auto_crop(seg)
% SEG_CROP()
%
% Automatic cropping of image segments.
% Input:
% seg - structure with informations about segment
% Output:
% seg_crop - cropped segment
%
% Author: Michal.Marczyk@polsl.pl

data = 255-seg.data;
[rows,cols] = find(data);
% Get the cropping parameters
topRow = min(rows(:));
bottomRow = max(rows(:));
leftColumn = min(cols(:));
rightColumn = max(cols(:));
% Extract a cropped image from the original.
data_crop = data(topRow:bottomRow, leftColumn:rightColumn);

seg_crop.x = seg.x(topRow:bottomRow);
seg_crop.y = seg.y(leftColumn:rightColumn);
seg_crop.data = 255 - data_crop;