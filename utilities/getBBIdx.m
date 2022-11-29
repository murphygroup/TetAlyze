% -------------------------------------------------------------------------
% Identifies coordinates for cortical BBs.
% -------------------------------------------------------------------------

function [cort_x, cort_y, cort_z, oaBBs] = getBBIdx(I, sensitivity, rejection_threshold, scale_xyz)
[height, width, zLength] = size(I);
[cell, cellArea] = identifyCell(I);
% if cellArea < 10000 % equivalent to 10000*(0.125^2) = 156.25 micrometer sq
%     [cell, ~] = identifyCell_triangle(I);
% end


[row, col] = find(cell == 1);

% boundaries enlarged by 0.125*40 = 5um
% image stack is cropped with the cell mask generated above (enlarged by
% 5um to account for irregularities in the mask)
enlg = 10;
minRow = max(1, min(row) - enlg);
minCol = max(1, min(col) - enlg);
maxRow = min(height, max(row) + enlg);
maxCol = min(width, max(col) + enlg);
I = I(minRow:maxRow, minCol:maxCol, :);
cell = cell(minRow:maxRow, minCol:maxCol);

[height, width, zLength] = size(I);


% convolved with a low radius (0.12um radius) Gaussian kernel
sigma = 1; 
I0 = imadjustn(I);
I_filtered = imgaussfilt3(I0, sigma);
T2 = 2.5*graythresh(nonzeros(I_filtered));
BW = imbinarize(I_filtered, T2);
SE = strel('sphere',1);
BW = imerode(BW, SE);
CC = bwconncomp(BW, 6);
S = regionprops(CC,'centroid', 'area');
OAarea = max([S.Area]);
L = labelmatrix(CC);
OAmask_ = ismember(L, find([S.Area] == OAarea));
OAmask = zeros(height, width, zLength);
OAmask(OAmask_)=1;
sliceViewer(OAmask);



sigma = 0.5; 

I = imadjustn(I);
% I_filtered = planeGaussianFilter(I, sigma);
I_filtered = imgaussfilt3(I, sigma);
I_filtered2 = imgaussfilt3(I, sigma);

% I_filtered(~BW) = 0;
% I_filtered2(~BW) = 0;
figure(1);
sliceViewer(I_filtered);

% I_filtered(I_filtered < meanT) = 0;

% T2 = 0.8*graythresh(nonzeros(I_filtered));
T2 = adaptthresh(I_filtered, sensitivity);
BW = imbinarize(I_filtered, T2);
BW =  bwareaopen(BW, 20);
% sliceViewer(BW);
SE = strel('sphere',1);
% SE = ones(3, 3, 3);
% SE3 = strel('disk', 1);
% SE2 = zeros(3, 3);
% SE2(2, 2) = 0.5;
% SE(1, :, :) = SE2;
% SE(2, :, :) = SE3.Neighborhood;
% SE(2, 1, 2) = 0.5;
% SE(2, 2, 1) = 0.5;
% SE(2, 3, 2) = 0.5;
% SE(2, 2, 3) = 0.5;
% SE(3, :, :) = SE2;
% sliceViewer(SE);
% BW = imerode(BW, SE);

% BW = imopen(BW, SE);
% BW = bwareaopen(BW,50);
% for i=1:zLength
%     A = BW(:,:,i);
%     A(~cell) =0;
%     BW(:,:,i) = A;
% end
% sliceViewer(BW);

I_filtered(~BW) = 0;
I_filtered2(~BW) = 0;

% BW = imregionalmax(I_filtered);
BW = activecontour(I_filtered2, BW);
BW = bwareaopen(BW,20);
figure(2);
sliceViewer(BW);

I_filtered(~BW) =0;
I_filtered2(~BW) =0;
CC = bwconncomp(BW);
s = regionprops(CC,'centroid', 'area');
L = labelmatrix(CC);

area = 100;
largeArea = ismember(L, find([s.Area] > area));
I_filtered2(~largeArea)=0;



BW = imregionalmax(I_filtered2);
idx = [s.Area] > area;

centers = cat(1, s.Centroid);
centers(idx,:) = [];
CC = bwconncomp(BW);
s2 = regionprops(CC,'centroid');
centers2 = cat(1,s2.Centroid);
% centers2 = [];
BW = vertcat(centers2, centers);

idx = [];


x = BW(:, 1)*0.108;
y = BW(:, 2)*0.108;
z = BW(:, 3)*0.2;
for i = 1; length(x)
    x_ = uint16(x(i)/0.108);
    y_ = uint16(y(i)/0.108);
    cube = cell(max(1, y_-1):min(y_+1, height), max(1, x_-1):min(x_+1, width));
    if sum(cube) == 0
        idx(end+1) = i;
    end
end
x(idx) = [];
y(idx) = [];
z(idx) = [];

% adaptive thresholding approach is used to separate noise from peaks
% corresponding to actual BB maxima
% To calculate the adaptive threshold, the average intensity and standard
% deviation for all peaks is calculated within a rolling average with a z
% range of 1.5um

% [y, x, z] = getPotentialBB_plane(I, BW);
% plot3(x,y,z,'ro')

% raw_points = [x(:),y(:),z(:)];

% distanceThreshold = 60; % unit is micron and is rouhly equal to the feret diameter constraint from identifyCell.m function
% To remove BB maxima within the interior of the cell, a 3D convex hull is
% generated from the BB maxima and all maxima that are greater than 2.25um
% from the surface of the convex hull are identified and deleted
shp = alphaShape(x, y, z);
shp.Alpha = 1.25*shp.Alpha;
[tri, xyz] = boundaryFacets(shp);
minDists = dist2AlphaShape(x, y, z, xyz, scale_xyz);
% minDists = dist2ConvexHull(x, y, z, k, scale_xyz);
% figure(9);
% plot(shp);
qualified = find(minDists <= rejection_threshold);

cort_x = x(qualified);
cort_y = y(qualified);
cort_z = z(qualified);


oaBBs = [];
for i = 1:length(cort_x)
    BB = [uint16(cort_x(i)/0.108), uint16(cort_y(i)/0.108), uint16(cort_z(i)/0.2)];
    cube = OAmask(max(1, BB(2)-1):min(BB(2)+1, height), max(1, BB(1)-1):min(BB(1)+1, width), max(1, BB(3)-1):min(BB(3)+1, zLength));
    if sum(cube) > 0
        oaBBs(end+1) = i;
    end

end    

end
% 
% [cort, oa] = distinguishBBClasses2(I, cort_x, cort_y, cort_z, OAmask);
% plot3(cort_x,cort_y,cort_z,'bo');


