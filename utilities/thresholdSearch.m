function [thresh, score] = thresholdSearch(I,zLength,scalar)

sigma = 0.5; 
% I_filtered = planeGaussianFilter(I, sigma);
I_filtered = imgaussfilt3(I, sigma);
I_filtered2 = imgaussfilt3(I, 2);
meanT = mean(I_filtered(:));
I_filtered2(I_filtered < meanT) = 0;
I_filtered(I_filtered < meanT) = 0;

thresh = scalar * graythresh(nonzeros(I_filtered));
BW = imbinarize(I_filtered, thresh);
I_filtered(~BW) =0;
I_filtered2(~BW) =0;
SE = strel('disk',1);
for i=1:zLength
   BW(:,:,i) = imerode(BW(:,:,i), SE);
end 
% BW = bwareaopen(BW,5);
I_filtered(~BW) = 0;
T = adaptthresh(I_filtered, 0.5);
BW = imbinarize(I_filtered, T);
I_filtered(~BW) =0;
I_filtered2(~BW) =0;

% BW = imregionalmax(I_filtered);
coefVar = std(nonzeros(I_filtered))/mean(nonzeros(I_filtered));

CC = bwconncomp(BW);
s = regionprops(CC,'centroid', 'area');

coefVar = std(nonzeros(I_filtered))/mean(nonzeros(I_filtered));
coefVar2 = std([s.Area])/mean([s.Area]);


score = coefVar2*coefVar;
end