function BW = restrict_mask(oriI)

% MIP created and mean pixel intensity is subtracted from MIP
[height, width, ~] = size(oriI);

oriMIP = sum(oriI, 3); % return maxes along dim 3 (depth)
totalMax = max(oriMIP,[],'all');
totalMin = min(oriMIP,[],'all');

% For grayscale images: For double arrays, values range from [0, 1]. For
% uint8, values range from [0, 255].
MIP = uint16(zeros(height, width));
MIP(:,:) = (oriMIP(:,:)-totalMin)/(totalMax-totalMin)*65536;

MIP_mean = mean(MIP,'all');
backgroundSubtractedMIP = MIP-(MIP_mean);

% Background subtracted MIP is convolved with large radius Gaussian kernel
% (1um radius) to get larger cellular scale features
backgroundSubtractedMIP = imgaussfilt(backgroundSubtractedMIP, 1); % standard deviation of 1

sigma = 1;
B = imgaussfilt(backgroundSubtractedMIP, sigma);
h = fspecial('log');
C = imfilter(B,h);
D = imgaussfilt(C, 10);


% Large features are separated by convolving with small radius Laplacian
% (0.12um radius) and further smoothed with another large radius Gaussian kernel
% h = fspecial('log'); % Laplacian of Gaussian filter
% C = imfilter(B, h);
% D = imgaussfilt(C, 25);


% resulting MIP is thresholded using Otsu's method and segmented objects
% are filtered based on their shape
level = graythresh(D); 
BW = imbinarize(D,level);
% imshow(BW);
end