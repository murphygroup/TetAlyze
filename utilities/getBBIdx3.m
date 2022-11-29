% -------------------------------------------------------------------------
% Identifies coordinates for cortical BBs.
% -------------------------------------------------------------------------

function [cort_x, cort_y, cort_z, antPole, postPole, order2OA, pot_oa, BW, I, Iori] = getBBIdx3(I, thresh_ratio, rejection_threshold, scale_xyz, resultPath, imageID, Iori)
% I = imgaussfilt3(I, 0.5);
res = [0.1806, 0.1806, 0.2002];
I = max(0, I - imgaussfilt3(I, 10));
I2 = fibermetric(I, 8);
[height, width, zLength] = size(I);
[cell, cellArea] = identifyCell(I2);

[row, col] = find(cell == 1);

% boundaries enlarged by 0.125*40 = 5um
% image stack is cropped with the cell mask generated above (enlarged by
% 5um to account for irregularities in the mask)
enlg = 3;
minRow = max(1, min(row) - enlg);
minCol = max(1, min(col) - enlg);
maxRow = min(height, max(row) + enlg);
maxCol = min(width, max(col) + enlg);
I = I(minRow:maxRow, minCol:maxCol, :);
Iori = Iori(minRow:maxRow, minCol:maxCol, :);
cell = cell(minRow:maxRow, minCol:maxCol);

OA_mask = getPotentialOARegion(I);

% imshow(cell);
SE = strel('disk',40);
cellS = imerode(cell,SE);
SE = strel('disk',10);
cellSS = imerode(cell,SE);
OA_mask = OA_mask .* cellSS;
% I = imresize3(I, 'Scale', [1, 1, 0.2/0.108]);
% Iori = imresize3(Iori, 'Scale', [1, 1, 0.2/0.108]);
[height, width, zLength] = size(I);

% 
% thresh = multithresh(I(:), 2);
% I(I < thresh(1)) = 0;



I2 = fibermetric(I, 4);
% sliceViewer(OAmask);
% I_filtered = imgaussfilt3(I2, 0.5);
OA_mask = getPotentialOARegion(I);
OA_mask = OA_mask .* cellSS;
% figure(11)
% imshow(OA_mask);
sigma = 0.5; 

I_filtered = I2;
I_filtered2 = I2; % imgaussfilt3(I, sigma);


thresh = 0.5* multithresh(I_filtered(:), 10);
BW = I_filtered > thresh(1+thresh_ratio-1);
if sum(BW(:))/(height*width*zLength) > 0.005
    BW = I_filtered > thresh(2+thresh_ratio-1);
end
if sum(BW(:))/(height*width*zLength) > 0.002
    BW = I_filtered > thresh(3+thresh_ratio-1);
end

% meanT = (1 + 0.2 * thresh_ratio) * mean(I_filtered(:));
% BW = I_filtered > meanT;


BW = activecontour(I_filtered, BW, 200, 'Chan-Vese', 'ContractionBias',-0.4);
% BW =  bwareaopen(BW, 20);

I_filtered(~BW) = 0;
I_filtered2(~BW) =0;



CC = bwconncomp(BW);
s = regionprops(CC,'centroid', 'area');
L = labelmatrix(CC);

area = 80;
largeArea = ismember(L, find([s.Area] > area));

CC = bwconncomp(largeArea);
% L = labelmatrix(CC);
% S = regionprops(CC,'centroid', 'area');

% SE = strel('sphere',1);
% largeArea = imerode(largeArea,SE);
I_filtered2(~largeArea)=0;


CC = bwconncomp(largeArea );
S2 = regionprops(CC,'centroid', 'area');
L = labelmatrix(CC);
areas = [S2.Area];
[areas, ~] = sort(areas, 'descend');
pot_oa2 = zeros(height, width, zLength);
pot_oa = zeros(height, width, zLength);
I5  = I;
for i=1:4
    [row,col] = find(sum(areas(i), 3)>0);
    a = zeros(height, width);
    for j=1:length(row)
        a(row(j), col(j)) = 1;
    end
    if areas(i) > 400 || sum(a .* cellS, 'all') > (0.9* sum(a, 'all'))
        oaArea = ismember(L, find([S2.Area] == areas(i))); 
        I_filtered2(oaArea)=0;
        pot_oa2 = pot_oa2 | oaArea;
        if i == 1
            I5(~pot_oa) = 0;
        end
%    
        figure(7);
        sliceViewer(oaArea);
    end
end

pot_oa_ori = imregionalmax(I.*pot_oa2);
CC = bwconncomp(pot_oa_ori);
s2 = regionprops(CC,'centroid');
pot_oa_ori = cat(1,s2.Centroid);

idx_ =  [];
a = prctile(I(largeArea),94);
for i=1:size(pot_oa_ori, 1)
    if I(pot_oa_ori(i, 2), pot_oa_ori(i, 1), pot_oa_ori(i, 3)) > a 
        idx_(end+1) = i;
    end
end
% pot_oa_ori_double = double(pot_oa_ori);



pot_oa =  zeros(length(idx_), 3);

pot_oa(:, 1) = pot_oa_ori(idx_, 1) * res(1);
pot_oa(:, 2) = pot_oa_ori(idx_, 2) * res(2);
pot_oa(:, 3) = pot_oa_ori(idx_, 3) * res(3);
pot_oa_idx = mean(pot_oa, 1);
if isempty(pot_oa)
    pot_oa_idx = mean(pot_oa_ori, 1);
    pot_oa_idx(1) = pot_oa_idx(1) * res(1);
    pot_oa_idx(2) = pot_oa_idx(1) * res(2);
    pot_oa_idx(3) = pot_oa_idx(1) * res(3);
end
[a, ~] = size(pot_oa);
pot_oa_ori2 = pot_oa_ori(idx_, :);

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

x = BW(:, 1)*res(1);
y = BW(:, 2)*res(2);
z = BW(:, 3)*res(3);

pot_oa_ori(idx_, :) = [];
x = [x' pot_oa_ori(:, 1)' * res(1)]';
y = [y' pot_oa_ori(:, 2)' * res(2)]';
z = [z' pot_oa_ori(:, 3)' * res(3)]';
BW = [BW' pot_oa_ori']';

% figure(7);
% sliceViewer(pot_oa_idx);
% 
% I_filtered2(oaArea)=0;


X = [x y z];
D = squareform(pdist(X));
for i=1:length(D)
    D(i, i) = D(i, i) + 100;
end
min_dists = min(D,[],2);
idx = find(min_dists > 3.75);
BW(idx, :) = []; 
x(idx, :) = []; 
y(idx, :) = []; 
z(idx, :) = []; 

shp = alphaShape(x, y, z);
shp.Alpha = 2*shp.Alpha;
[tri, xyz] = boundaryFacets(shp);
minDists = dist2AlphaShape(x, y, z, xyz, scale_xyz);
% minDists = dist2ConvexHull(x, y, z, k, scale_xyz);
% figure(9);
% plot(shp);
qualified = find(minDists <= rejection_threshold);

BW = BW(qualified, :);


cort_x = BW(:, 1)*res(1);
cort_y = BW(:, 2)*res(2);
cort_z = BW(:, 3)*res(3);



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

% rescue BBs mistakenly identified as OA BBs
[idx1, ~] = find(pot_oa(:, 3) > max(pot_oa(:, 3))*0.95);
cort_x = [cort_x' pot_oa(idx1, 1)']';
cort_y = [cort_y' pot_oa(idx1, 2)']';
cort_z = [cort_z' pot_oa(idx1, 3)']';

pot_oa(idx1, :) = [];
% 
BW = [BW' pot_oa_ori2(idx1, :)']';


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

for i=1:zLength
    idx = find(z == i);
    Xs = BW(idx, 1);
    Ys = BW(idx, 2);
    if ~isempty(idx)
        if i < startFrame
            startFrame = i;
        end
        if i > endFrame
            endFrame = i;
        end
        img = I(:, :, i)*10;
        img = cat(3, img, img, img);
        imgStack = cat(4, imgStack, img);
        for j =1:length(idx)
            img = insertMarker(img,[Xs(j) Ys(j)],'circle', 'size',4);
        end

        imgStack = cat(4, imgStack, img);
    end
end
p = 0;
for i=startFrame+1:zLength
    idx = find(z == i);
    Xs = BW(idx, 1);
    Ys = BW(idx, 2);
    if ~isempty(idx)
%         img = imgStack(:, :, :, 2*p); 
        p = p+1;
        for j =1:length(idx)
            imgStack(:, :, :, 2*p) = insertMarker(imgStack(:, :, :, 2*p),[Xs(j) Ys(j)],'star', 'size',3, 'Color', 'red');
        end
    end
end
p = 1;
for i=1:endFrame-1
    idx = find(z == i);
    Xs = BW(idx, 1);
    Ys = BW(idx, 2);
    if ~isempty(idx)
%         img = imgStack(:, :, :, 2*p); 
        p = p+1;
        for j =1:length(idx)
            imgStack(:, :, :, 2*p) = insertMarker(imgStack(:, :, :, 2*p),[Xs(j) Ys(j)],'star', 'size',3, 'Color', 'cyan');
        end
    end
end
imgStack = permute(imgStack, [1 2 4 3]);
% a = size(imgStack)
f = figure(1);
sliceViewer(imgStack);
save(join([resultPath, imageID, 'BB Identification.mat']), 'imgStack');

if isempty(lowestOA)
    lowestOA = pot_oa_idx;
end
dist2OA = vecnorm([cort_x cort_y cort_z] - lowestOA(1, :),2,2);
[~, order2OA] = sort(dist2OA,'ascend');
end