function [x,y,data] = auto_crop2(x,y,data)
% SEG_CROP()
%
% Automatic cropping of image segments.
% Input:
% seg - structure with informations about segment
% Output:
% seg_crop - cropped segment
%
% Author: Michal.Marczyk@polsl.pl

data = data;
[rows,cols] = find(data);

% Get the cropping parameters
topRow = min(rows(:));
bottomRow = max(rows(:));
leftColumn = min(cols(:));
rightColumn = max(cols(:));

% Extract a cropped image from the original.
data = data(topRow:bottomRow, leftColumn:rightColumn);
x = x(topRow:bottomRow);
y = y(leftColumn:rightColumn);
