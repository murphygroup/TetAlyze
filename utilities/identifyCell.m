% -------------------------------------------------------------------------
%% [Mehul] 18-Feb-2020 - Some general notes to improve work flow understanding

    % (Processed image data in 3-D array) 
    %   --> (Apply min-max normalization)
    %   --> (Substract Mean to remove noise)
    %   --> Apply Gaussian Filter with sigma(standard deviation) = 1
    %   --> Apply Laplacian of Gaussian Filter
    %   --> Apply Gaussian Filter with sigma(standard deviation) = 8 for
    %       extra smoothening
    %   --> Apply Otsu's global thresholding technique to distinguish
    %       background and foreground
    %   --> Binarizing the image using the Otsu's threshold obtained from
    %       above
    %   Some constraints which are hard coded are:
    %   --> area of each cell <=128000
    %   --> feret diameter < 60 micron
    %   --> circularity >= 0.85
                        
                        
%   [Ben] (adapted from Jingyi's code by Ben) 01/02/18
% Performs the cell detection outlined in pg 16 of Jinyi's paper.
% Extracting MIP, applying several filters, binarizing image, searching for
% largest object in resultant image satisfying several conditions.
% Returns mask of detected cell, together with area of mask.
% -------------------------------------------------------------------------

%%
function [cell,cellArea] = identifyCell(oriI)

% MIP created and mean pixel intensity is subtracted from MIP
[height, width, ~] = size(oriI);

oriMIP = max(oriI, [], 3); % return maxes along dim 3 (depth)
totalMax = max(oriMIP,[],'all');
totalMin = min(oriMIP,[],'all');

% For grayscale images: For double arrays, values range from [0, 1]. For
% uint8, values range from [0, 255].
MIP = uint16(zeros(height, width));
MIP(:,:) = (oriMIP(:,:)-totalMin)/(totalMax-totalMin)*65536;

MIP_mean = 0.1*mean(MIP,'all');
backgroundSubtractedMIP = MIP-(MIP_mean);

% Background subtracted MIP is convolved with large radius Gaussian kernel
% (1um radius) to get larger cellular scale features
B = imgaussfilt(backgroundSubtractedMIP, 1); % standard deviation of 1
% BW = B > mean(B(:));
% B(~BW) = 0;

% B = imgaussfilt(backgroundSubtractedMIP, sigma, 'FilterSize', 1/0.125);
% h = fspecial('log', 0.12/0.125, sigma); 
% D = imgaussfilt(C, sigma, 'FilterSize', 1/0.125);


% Large features are separated by convolving with small radius Laplacian
% (0.12um radius) and further smoothed with another large radius Gaussian kernel
h = fspecial('log'); % Laplacian of Gaussian filter
C = imfilter(B, h);
D = imgaussfilt(C, 2.5);


% resulting MIP is thresholded using Otsu's method and segmented objects
% are filtered based on their shape
level = graythresh(D); 
BW = imbinarize(D, level);
se = strel('disk', 8);
BW = imdilate(BW,se);
BW = ~bwareaopen(~BW, 5000);
% imshow(BW);
L = bwlabel(BW); % label connected components
areaSize = zeros( max(L,[],'all'), 1);

for label_num = 1:max(L,[],'all')
    [y, x] = find(L == label_num);
    try
        [~, areaSize(label_num)] = convhull(x,y);
        % areaSize(label_num)=polyarea(x(K), y(K)); (old code for new
        % images. The above updated one directly calculates the area.
        % Area represented by 1 pixel seems to be (0.125^2) 
        % objects > 2000um^2 in 2D area
        
%         if areaSize(label_num) > 4000/(0.108^2)
%             areaSize(label_num) = 0;
%         end
% 
%         feret diameter > 60um
        currImg = zeros(height, width);
        for i = 1:length(x)
            currImg(y(i), x(i)) = 1;
        end
        
        currImg = imfill(currImg);
        [fd, ~] = imFeretDiameter(currImg);
%         if fd > 100/0.108
%             areaSize(label_num) = 0;
%         end
        
%         circularity < 0.75
        stats = regionprops(currImg, 'Perimeter', 'Area');   
  
%         circularity = 4*pi*stats.Area/stats.Perimeter^2;
%         if circularity < 0.75 % perfect circle has circularity == 1
%             areaSize(label_num) = 1;
%         end 
%     catch
%         areaSize(label_num) = 0;
    end
end

% largest remaining object that is retained is the cell of interest
[cellArea, label] = max(areaSize);
% label = find(areaSize == cellArea, 1);

[y, x] = find(L == label);
% k1 = convhull(x,y);
% plot(x(k1),y(k1))
% hold on
cell = zeros(height, width);

for i=1:length(x)
    cell(y(i), x(i)) = 1;
end

cell = imfill(cell);
se = strel('disk',7);
cell = imopen(cell,se);
% hold on;
% imshow(cell)
end