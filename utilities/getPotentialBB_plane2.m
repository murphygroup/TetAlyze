% -------------------------------------------------------------------------
% [Ben] 12/07/17 
% Takes in cropped image stack, corresponding binary array indicating
% location of local maxima, and parameter used in calculation of adaptive
% threshold.
% Further filters local maxima, and returns coordinates of remaining
% potential BBs.
% -------------------------------------------------------------------------

function [y, x, z] = getPotentialBB_plane2(Img, BW, th)
[m, n, l] = size(Img);

%Img_vect = reshape(Img,m*n*l,1);
%Img_vect = sort(Img_vect);

thresholds = zeros(l, 1);
%size(Img_vect)
%my_thresh = Img_vect(round(m*n*l*0.997),1)
target_object_count = 400
decrease = 0.8

brightness_decrease_factor = zeros(l)
for i = 1:l
    brightness_decrease_factor(i) = 1.0-(i-1)*(1.0-decrease)/i

for i = 1:l
    % 'found local maxima within a box with size approximating the size of
    % basal body 0.25x0.25x0.6 (um)^3'
    % i.e. each plane taken to have thickness of 3~5 pixels
    if i == 1
        currPlanes = Img(:, :, i:i+2);
        currBWPlanes = BW(:, :, i:i+2);
        % row, col, depth coordinates of identified regional maxima
        % note that 'find(X)' returns linear indices
        [row, col, depth] =  ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
    elseif i == 2
        currPlanes = Img(:, :, i-1:i+2);
        currBWPlanes = BW(:, :, i-1:i+2);
        [row, col, depth] = ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
    elseif i == l-1
        currPlanes = Img(:, :, i-2:i+1);
        currBWPlanes = BW(:, :, i-2:i+1);
        [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
    elseif i==l
        currPlanes=Img(:, :, i-2:i);
        currBWPlanes=BW(:, :, i-2:i);
        [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
    else
        currPlanes=Img(:, :, i-2:i+2);
        currBWPlanes=BW(:, :, i-2:i+2);
        [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
    end
    
    currPlaneIntensity = zeros(length(row), 1);
    for pt = 1:length(row)
        currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
    end
    
    %currMean = mean(currPlaneIntensity)
    currMean = mean(mean(mean(currPlanes)))
    mean(mean(mean(currBWPlanes)))
    
    currStd = std(currPlaneIntensity)
    
    currCv = currStd/currMean
    
    % pg. 17 of Jingyi's paper: 'intensity threshold for given plane ...'
    currThreshold = (th - th*currCv)*currStd + currMean % given in Chad's paper
    %((6?(Icv×6))×Is.d.)+Iavg
    currThreshold = (th-(currCv*th))*currStd + currMean
    
    % thresholds(i) = currThreshold;
    %prctile(Img,1,'all')
    %thresholds(i) = prctile(Img,99,'all')
    %my_thresh = graythresh(currPlanes)
    %thresholds(i) = prctile(currPlaneIntensity,my_thresh,'all')
    my_max = max(max(max(Img)));
    %size(my_thresh)
    %size(my_max)
    %thresholds(i) = currMean + 2 * currStd
    mean(currPlanes);
    thresholds(i) = currThreshold%mean(mean(mean(currPlanes)))*50;
    %thresholds(i) = mean(mean(mean(currPlanes)))*10;
    %thresholds(i) = max(max(max(currPlanes))) * mean(mean(mean(currPlanes)))
    %thresholds(i,1) = currCv
    %thresholds(i,2) = currMean;
    %thresholds(i,3) = currStd
    %thresholds(i,4) = currThreshold
    
end
thresholds
plot(thresholds)
mtx = zeros(size(BW));

sum(sum(sum(BW)))


for i=1:m
    for j=1:n
        for k=1:l
            % adaptive threshold used to filter local maxima
            if BW(i, j, k) == 1 && Img(i, j, k) >= thresholds(k)
                mtx(i, j, k) = Img(i, j, k);
            end
        end
    end
end

[y, x, z] = ind2sub(size(mtx), find(mtx > 0));

size(y)
end
