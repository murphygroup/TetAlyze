% -------------------------------------------------------------------------
%% [Mehul] 5-March-2020
% This function is similar to identifyCell function. It is just applying
% the gaussian filters along with 'log' and unsharp filter in a loop for
% better detection of basal bodies in a football shape.
%%
function BW = calibrate_detect(oriI, z)

[p,q,r] = size(oriI);
% max_tol = 0.2;
% H_val = 1 - max_tol;
% fmaxtol = @(x) x-(max_tol*x);
npc = 4;

switch(z)
    case 1
%         z_init = z;
%         z_fin = z_init+npc;
        currPlanes = max(oriI(:, :, z:z+2),[],3);
        cp_max = max(currPlanes,[],'all');
        cp_min = min(currPlanes,[],'all');
        cpMIP = double(zeros(p, q));
        cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);
        
        
    case 2
%         z_init = z;
%         z_fin = z_init+npc;
        currPlanes = max(oriI(:, :, z-1:z+2),[],3);
        cp_max = max(currPlanes,[],'all');
        cp_min = min(currPlanes,[],'all');
        cpMIP = double(zeros(p, q));
        cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);
        
    case r-1
        z_init = z-npc+1;
        z_fin = z_init+npc;
        currPlanes = max(oriI(:, :, z-2:z+1),[],3);
        cp_max = max(currPlanes,[],'all');
        cp_min = min(currPlanes,[],'all');
        cpMIP = double(zeros(p, q));
        cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);
        
    case r
%         z_init = z-npc;
%         z_fin = z_init+npc;
        currPlanes = max(oriI(:, :, z-2:z),[],3);
        cp_max = max(currPlanes,[],'all');
        cp_min = min(currPlanes,[],'all');
        cpMIP = double(zeros(p, q));
        cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);
        
%     case {floor(r/2),floor(r/2)+1}
%         z_init = z-1;
%         z_fin = z_init+npc;
%         currPlanes = max(oriI(:, :, z_init:z_fin),[],3);
%         cp_max = max(currPlanes,[],'all');
%         cp_min = min(currPlanes,[],'all');
%         cpMIP = double(zeros(p, q));
%         cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);
%         
%     case {r-2,r-3,r-4,r-5}
%         z_init = z-npc;
%         z_fin = z_init+npc;
%         currPlanes = max(oriI(:, :, z_init:z_fin),[],3);
%         cp_max = max(currPlanes,[],'all');
%         cp_min = min(currPlanes,[],'all');
%         cpMIP = double(zeros(p, q));
%         cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);

    otherwise
%         z_init = z;
%         z_fin = z_init+npc;
        currPlanes = max(oriI(:, :, z-2:z+2),[],3);
        cp_max = max(currPlanes,[],'all');
        cp_min = min(currPlanes,[],'all');
        cpMIP = double(zeros(p, q));
        cpMIP(:,:) = (currPlanes(:,:)-cp_min)/(cp_max-cp_min);
end

cpMIP = cpMIP - mean(cpMIP,'all');
% cpMIP = imsharpen(cpMIP);
level = graythresh(cpMIP);
BW = imbinarize(cpMIP, level);
% size(BW)
% imshow(BW);

end


% for i=1:ceil(p1)
%     B = imgaussfilt(B_temp, p2 ); 
%     h = fspecial('motion',p3,p4); 
%     B_temp = imfilter(B, h);
%     h1 = fspecial('log'); % Laplacian of Gaussian filter
%     B_temp = imfilter(B_temp,h1);
%     D = imgaussfilt(B_temp,p5);
% end

























% % aall = [13.0739    2.2054    3.2566   25.0706   23.3739];
% % p1 = aall(1);
% % p2 = aall(2);
% % p3 = aall(3);
% % p4 = aall(4);
% % p5 = aall(5);
% 
% % MIP created and mean pixel intensity is subtracted from MIP
% [height, width, depth] = size(oriI);
% 
% oriMIP = max(oriI, [], 3); % return maxes along dim 3 (depth)
% totalMax = max(oriMIP,[],'all');
% totalMin = min(oriMIP,[],'all');
% 
% % For grayscale images: For double arrays, values range from [0, 1]. For
% % uint8, values range from [0, 255].
% MIP = uint16(zeros(height, width));
% MIP(:,:) = (oriMIP(:,:)-totalMin)/(totalMax-totalMin)*65536;
% 
% MIP_mean = mean(MIP,'all');
% backgroundSubtractedMIP = MIP-(MIP_mean);
% 
% % Background subtracted MIP is convolved with large radius Gaussian kernel
% % (1um radius) to get larger cellular scale features
% B = imgaussfilt(backgroundSubtractedMIP, 1); % standard deviation of 1
% 
% h = fspecial('log'); % Laplacian of Gaussian filter
% C = imfilter(B, h);
% D = imgaussfilt(C, 8);
% 
% 
% % resulting MIP is thresholded using Otsu's method and segmented objects
% % are filtered based on their shape
% level = graythresh(D); 
% BW = imbinarize(D, level);
% 
% L = bwlabel(BW); % label connected components
% areaSize = zeros( max(L,[],'all'), 1);
% 
% for label_num = 1:max(L,[],'all')
%     [y, x] = find(L == label_num);
%     try
%         [~, areaSize(label_num)] = convhull(x,y);
%         % areaSize(label_num)=polyarea(x(K), y(K)); (old code for new
%         % images. The above updated one directly calculates the area.
%         % Area represented by 1 pixel seems to be (0.125^2) 
%         % objects > 2000um^2 in 2D area
%         
%         if areaSize(label_num) > 2000/(0.125^2)
%             areaSize(label_num) = 0;
%         end
% 
%         % feret diameter > 60um
%         currImg = zeros(height, width);
%         for i = 1:length(x)
%             currImg(y(i), x(i)) = 1;
%         end
%         
%         currImg = imfill(currImg);
%         [fd, ~] = imFeretDiameter(currImg);
%         if fd > 60/0.125
%             areaSize(label_num) = 0;
%         end
%         
%         % circularity < 0.85
%         stats = regionprops(currImg, 'Perimeter', 'Area');   
%   
%         circularity = 4*pi*stats.Area/stats.Perimeter^2;
%         if circularity < 0.85 % perfect circle has circularity == 1
%             areaSize(label_num) = 0;
%         end 
%     catch
%         areaSize(label_num) = 0;
%     end
% end
% 
% % largest remaining object that is retained is the cell of interest
% [cellArea, label] = max(areaSize);
% % label = find(areaSize == cellArea, 1);
% 
% [y, x] = find(L == label);
% % k1 = convhull(x,y);
% % plot(x(k1),y(k1))
% % hold on
% cell = zeros(height, width);
% 
% for i=1:length(x)
%     cell(y(i), x(i)) = 1;
% end
% 
% cell = imfill(cell);
% % hold on;
% % imshow(cell)
% end

