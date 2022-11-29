% -------------------------------------------------------------------------
% Identifies coordinates for cortical BBs.
% -------------------------------------------------------------------------

function [cort_x, cort_y, cort_z, antPole, postPole, order2OA, pot_oa, BW, I] = getBBIdx_4correction(I, rejection_threshold, scale_xyz, resultPath, imageID, Iori)
I2 = fibermetric(I, 8);
[height, width, zLength] = size(I);
[cell, cellArea] = identifyCell(I2);

[row, col] = find(cell == 1);

% boundaries enlarged by 0.125*40 = 5um
% image stack is cropped with the cell mask generated above (enlarged by
% 5um to account for irregularities in the mask)
enlg = 20;
minRow = max(1, min(row) - enlg);
minCol = max(1, min(col) - enlg);
maxRow = min(height, max(row) + enlg);
maxCol = min(width, max(col) + enlg);
I = I(minRow:maxRow, minCol:maxCol, :);
cell = cell(minRow:maxRow, minCol:maxCol);

[height, width, zLength] = size(I);
% imshow(cell);


I2 = fibermetric(I, 8);
% sliceViewer(OAmask);
% I_filtered = imgaussfilt3(I2, 0.5);
OA_mask = getPotentialOARegion(I);
% figure(11)
% imshow(OA_mask);
sigma = 1; 

I_filtered = I2;
I_filtered2 = imgaussfilt3(I, sigma);


meanT = 0.5*mean(I_filtered(:));
BW = I_filtered > meanT;

BW = activecontour(I_filtered, BW);
BW =  bwareaopen(BW, 3);

I_filtered(~BW) =0;
I_filtered2(~BW) =0;


CC = bwconncomp(BW);
s = regionprops(CC,'centroid', 'area');
L = labelmatrix(CC);

area = 80;
largeArea = ismember(L, find([s.Area] > area));

CC = bwconncomp(largeArea);
% L = labelmatrix(CC);
% S = regionprops(CC,'centroid', 'area');

SE = strel('sphere',1);
largeArea = imerode(largeArea,SE);
I_filtered2(~largeArea)=0;


CC = bwconncomp(largeArea );
S2 = regionprops(CC,'centroid', 'area');
L = labelmatrix(CC);
areas = [S2.Area];
[areas, ~] = sort(areas, 'descend');
pot_oa = zeros(height, width, zLength);
I5  = I;
for i=1:2
    if areas(i) > 75 
        oaArea = ismember(L, find([S2.Area] == areas(i)));
        I_filtered2(oaArea)=0;
%         pot_oa = pot_oa | oaArea;
        if i == 1
            pot_oa2 = oaArea;
            I5(~pot_oa2) = 0;
        end
   
%         figure(7);
%         sliceViewer(oaArea);
    end
end
% figure(9);
% sliceViewer(I5);
pot_oa = imregionalmax(I5);
CC = bwconncomp(pot_oa);
s2 = regionprops(CC,'centroid');
pot_oa = cat(1,s2.Centroid);

pot_oa(:, 1) = pot_oa(:, 1) * 0.108;
pot_oa(:, 2) = pot_oa(:, 2) * 0.108;
pot_oa(:, 3) = pot_oa(:, 3) * 0.2;
pot_oa_idx = mean(pot_oa, 1);
[a, ~] = size(pot_oa);

% figure(7);
% sliceViewer(pot_oa_idx);
% 
% I_filtered2(oaArea)=0;

BW = imregionalmax(I_filtered2);
idx = [s.Area] > area;

centers = cat(1, s.Centroid);
centers(idx,:) = [];
CC = bwconncomp(BW);
s2 = regionprops(CC,'centroid');
centers2 = cat(1,s2.Centroid);
% centers2 = [];

BW = vertcat(centers2, centers);
% 
idx = [];
for i=1:length(BW)
    x_ = cast(BW(i, 1), 'uint32');
    y_ = cast(BW(i, 2), 'uint32');
    if cell(y_, x_) == 0
        idx(end+1) = i;
    end
end
BW(idx, :) = []; 

x = BW(:, 1)*0.108;
y = BW(:, 2)*0.108;
z = BW(:, 3)*0.2;


shp = alphaShape(x, y, z);
shp.Alpha = 1.25*shp.Alpha;
[tri, xyz] = boundaryFacets(shp);
minDists = dist2AlphaShape(x, y, z, xyz, scale_xyz);
% minDists = dist2ConvexHull(x, y, z, k, scale_xyz);
% figure(9);
% plot(shp);
qualified = find(minDists <= rejection_threshold);

BW = BW(qualified, :);
[~, oaIdx] =  distinguishBBClasses2(I, BW(:, 1), BW(:, 2), BW(:, 3), OA_mask);

BW(oaIdx, :) = [];

cort_x = BW(:, 1)*0.108;
cort_y = BW(:, 2)*0.108;
cort_z = BW(:, 3)*0.2;



[pole1, pole2] = findPoles(cort_x, cort_y, cort_z,scale_xyz);
% oa_x = cort_x(pot_oa_idx);
% oa_y = cort_y(pot_oa_idx);
% oa_z = cort_z(pot_oa_idx);
% OA = round([mean(oa_x), mean(oa_y), mean(oa_z)]);
d1 = distance_pts(pole1, pot_oa_idx, [1,1,1]);
d2 = distance_pts(pole2, pot_oa_idx, [1,1,1]);
% anterior pole is closer to OA region
if d1 > d2
    antPole = pole2;
    postPole = pole1;
else
    antPole = pole1;
    postPole = pole2;
end
cort_x(end+1:end+a) = pot_oa(:, 1);
cort_y(end+1:end+a) = pot_oa(:, 2);
cort_z(end+1:end+a) = pot_oa(:, 3);

[cort_x, cort_y, cort_z] = newCoorWithoutOA(cort_x, cort_y, cort_z, antPole, postPole);
antPole = [0 0 norm(pole1 - pole2)];
postPole = [0 0 0];

pot_oa(:, 1) = cort_x(end-a+1:end);
pot_oa(:, 2) = cort_y(end-a+1:end);
pot_oa(:, 3) = cort_z(end-a+1:end);
idx = find(pot_oa(:, 3) == min(pot_oa(:, 3)));

lowestOA = [pot_oa(idx, 1), pot_oa(idx, 2), pot_oa(idx, 3)];
cort_x(end-a+1:end) = [];
cort_y(end-a+1:end) = [];
cort_z(end-a+1:end) = [];

% if ~isempty(oaIdx)
%     BW(oaIdx, :) = [];
%     cort_x(oaIdx)=[];
%     cort_y(oaIdx)=[];
%     cort_z(oaIdx)=[];
% end



qualified = ones(length(cort_x),1);
[idx1, ~] = find(cort_z > max(cort_z) - 13);
[idx2, ~] = find(cort_z < max(cort_z) - 3);
idx = intersect(idx1, idx2);
shp1 = alphaShape(0.9*cort_x(idx), 0.9*cort_y(idx), cort_z(idx));
if shp1.Alpha ~= inf
    BBs = [0.5*cort_x(idx), 0.5*cort_y(idx), cort_z(idx)];
    shp1 = alphaShape(BBs(:, 1), BBs(:, 2), BBs(:, 3));
    shp1.Alpha = inf;
    tf = inShape(shp1, cort_x(idx), cort_y(idx), cort_z(idx));
%     figure(6);
%     plot(shp1);

else
    tf = ones(length(idx), 1);
end
for j = 1:length(idx)
    if tf(j) ~= 0 
        qualified(idx(j)) = 0;
    end
end
cort_x = cort_x(find(qualified));
cort_y = cort_y(find(qualified));
cort_z = cort_z(find(qualified));

BW = BW(find(qualified), :);

imgStack = [];
startFrame = 1000000;
endFrame = -1;
z = cast(BW(:, 3), 'uint32');

dist2OA = vecnorm([cort_x cort_y cort_z] - lowestOA(1, :),2,2);
[~, order2OA] = sort(dist2OA,'ascend');
end