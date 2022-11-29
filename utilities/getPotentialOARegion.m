% -------------------------------------------------------------------------
% [Ben] 3/16/18 (adapted from Jingyi's getOA.m)
% Gets the coordinates describing the potential OA region in pg 17 of
% Jingyi's paper. These coordinates are stored in the form of a binary
% matrix, and will be subsequently be used in the differentation of
% cortical and oral apparatus BBs.
% -------------------------------------------------------------------------


function OA_mask = getPotentialOARegion(I)
% 'summed all intensities at each position for all images within image
% stack' pg 17.

SIP = I; % sum(I, 3); % summed intensity projection
minInt = min(SIP,[],'all');
maxInt = max(SIP,[],'all');
ratio = 0.35;

% 'set a threshold for intensity and a threshold for area' pg 17.
th_int = minInt + ratio*(maxInt - minInt);
% th_int = graythresh(SIP);
th_area = 50; 

BW = (SIP > th_int)*1;

[L, n] = bwlabeln(BW); % n is # of connected components
nums = zeros(n, 1);
% figure(13);
% imshow(BW);

for label=1:n
    % nums stores the 'sizes' of the various connected components
    nums(label) = length(find(L == label));
end
areas = maxk(nums, 3);
for label=1:n
    if nums(label) < th_area
        BW(find(L == label)) = 0;
    end
end
se = strel('disk', 5);
BW = imdilate(BW,se);
% r stores the indices (relative to nums) of the connected components that
% are 'large' enough
% r = find(nums > th_area);

[L, n] = bwlabeln(BW); % n is # of connected components
nums = zeros(n, 1);
for label=1:n
    % nums stores the 'sizes' of the various connected components
    nums(label) = length(find(L == label));
end

if n > 3
    [~, idx] = maxk(nums, 4); 
    BW(find(L ~= idx(1) & L ~= idx(2) & L ~= idx(3) & L ~= idx(4))) = 0;
    OA_mask = BW;
elseif n > 2
    [~, idx] = maxk(nums, 3); 
    BW(find(L ~= idx(1) & L ~= idx(2) & L ~= idx(3) )) = 0;
    OA_mask = BW;
elseif n > 1
    [~, idx] = maxk(nums, 2); 
    BW(find(L ~= idx(1) & L ~= idx(2))) = 0;
    OA_mask = BW;
elseif n > 0
    [~, idx] = max(nums); 
    BW(find(L ~= idx)) = 0;

    % imshow(BW);
    OA_mask = BW;
else
    [a, b, c] = size(SIP);
    OA_mask = zeros(a, b, c);
end

% for i=1:length(r)
%     [row, col] = find(L == r(i));
%     % don't correct the indices for cropping here
%     for j=1:length(row)
%         OA_mask(row(j), col(j)) = 1;
%     end
% end
end

