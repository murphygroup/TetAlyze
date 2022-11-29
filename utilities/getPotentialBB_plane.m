% -------------------------------------------------------------------------
%% Mehul - 13 March 2020 - Replaced with modified code with a modified version of identifyCell
    % implemented in this one. The intensity thresholds are changed to
    % intensities having 1/8th intensity of the maximum intensity in that
    % particular plane.

% [Ben] 12/07/17 
% Takes in cropped image stack, corresponding binary array indicating
% location of local maxima, and parameter used in calculation of adaptive
% threshold.
% Further filters local maxima, and returns coordinates of remaining
% potential BBs.
% -------------------------------------------------------------------------

function [y, x, z] = getPotentialBB_plane(Img, BW)
[p, q, r] = size(Img);
% th = 1;
% Img = bw01(Img);
% max(Img,[],'all')
p1 = 1.3;
p2 = 2.7;

% Centroids = [3.1164   46.8690;9.9307  207.8219; 7.2096  127.3564];
allmin = zeros(1, r);
% sk_Img = skewness(Img,1,[1,2]);
% avg_sk = mean(sk_Img(:));
% ku_Img = kurtosis(Img,1,[1,2]);
% avg_ku = mean(ku_Img(:));
% 
% % thfunc = @(x) mean(x);
% 
% if avg_sk > 5.5 && avg_ku > 100
%     thfunc = @(x) mean(x)*p1;
% else
%     thfunc = @(x) mean(x)*p2;
% end




for i = 1:r
    sk_Img = skewness(Img(:,:,i),1,[1,2]);
    ku_Img = kurtosis(Img(:,:,i),1,[1,2]);
    med_Img = median(Img(:,:,i),[1,2]);
    std_Img = std(Img(:,:,i),0,[1,2]);
    
    p_val_fit = [-5.41935065019943,0.124000181372767,-0.0621734206226516,-0.148696503642740,0.0868141547505692,0.521167418810021,4.58288128339685,0.0121263833819451,0.00457541599613812,4.61437989184478];
    p_val = p_val_fit(1)*(ku_Img/sk_Img)^p_val_fit(2)...
        + p_val_fit(3)*(sin(p_val_fit(4)*ku_Img)) + p_val_fit(5)*cos(p_val_fit(6)*ku_Img) ...
        + p_val_fit(7)*med_Img^p_val_fit(8)/std_Img^p_val_fit(9)  + p_val_fit(10);

    thfunc = @(x) mean(x)*p_val;
    
     % 'found local maxima within a box with size approximating the size of
    % basal body 0.25x0.25x0.6 (um)^3'
    % i.e. each plane taken to have thickness of 3~5 pixels
    switch(i)
    case 1
        currPlanes = Img(:, :, i:i+2);
        currBWPlanes = BW(:, :, i:i+2);
        [row, col, depth] =  ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
        currPlaneIntensity = zeros(length(row), 1);
        for pt = 1:length(row)
            currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
        end
%         allmin_temp = sort(currPlaneIntensity,'descend');
        allmin(i) = thfunc(currPlaneIntensity);%allmin_temp(thr_pts(i));%
        
               
    case 2
        currPlanes = Img(:, :, i-1:i+2);
        currBWPlanes = BW(:, :, i-1:i+2);
        [row, col, depth] = ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
        currPlaneIntensity = zeros(length(row), 1);
        for pt = 1:length(row)
            currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
        end
%         allmin_temp = sort(currPlaneIntensity,'descend');
        allmin(i) = thfunc(currPlaneIntensity);%allmin_temp(thr_pts(i));%
        
    case r-1
        currPlanes = Img(:, :, i-2:i+1);
        currBWPlanes = BW(:, :, i-2:i+1);
        [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
        currPlaneIntensity = zeros(length(row), 1);
        for pt = 1:length(row)
            currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
        end
%         allmin_temp = sort(currPlaneIntensity,'descend');
        allmin(i) = thfunc(currPlaneIntensity);%allmin_temp(thr_pts(i));%
        
    case r
        currPlanes=Img(:, :, i-2:i);
        currBWPlanes=BW(:, :, i-2:i);
        [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
        currPlaneIntensity = zeros(length(row), 1);
        for pt = 1:length(row)
            currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
        end
%         allmin_temp = sort(currPlaneIntensity,'descend');
        allmin(i) = thfunc(currPlaneIntensity);%allmin_temp(thr_pts(i));%
        
    otherwise
        currPlanes=Img(:, :, i-2:i+2);
        currBWPlanes=BW(:, :, i-2:i+2);
        [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
        currPlaneIntensity = zeros(length(row), 1);
        for pt = 1:length(row)
            currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
        end
%         allmin_temp = sort(currPlaneIntensity,'descend');
        allmin(i) = thfunc(currPlaneIntensity);%allmin_temp(thr_pts(i));%
    end

end
% figure(5)
% plot(allmin,'ro')

yfit_data = allmin;
% xfit_data = (1:r);
% 
% 
thresholds = yfit_data; %threshfit(xfit_data,yfit_data); %yfit_data;
% % hold on
% % plot(thresholds,'g.')
% % hold off


mtx = zeros(size(BW));

[cell1,~] = identifyCell(Img);



for k=1:r
    cell2 = calibrate_detect(Img,k);
%     disp(strcat('At plane : ',num2str(k)))
    for j=1:q
        for i=1:p
            % adaptive threshold used to filter local maxima
            if BW(i, j, k) == 1 && Img(i, j, k) >= thresholds(k)
                if cell2(i,j)==1 && cell1(i,j) ==1
                    mtx(i, j, k) = Img(i, j, k);
                end
            end
        end
    end
end



[y, x, z] = ind2sub(size(mtx), find(mtx > 0));
end






















%     case {floor(r/2),floor(r/2)+1}
%         currPlanes=Img(:, :, i-1:i+1);
%         currBWPlanes=BW(:, :, i-1:i+1);
%         [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%         currPlaneIntensity = zeros(length(row), 1);
%         for pt = 1:length(row)
%             currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
%         end
% %         allmin_temp = sort(currPlaneIntensity,'descend');
%         allmin(i) = mean(currPlaneIntensity);%allmin_temp(thr_pts(i));%


% a1 =       31.84;
% b1 =     0.06424;
% c1 =     0.06522;
% a2 =         6.4;
% b2 =      0.2518;
% c2 =       1.721;
% thr_pts = ones(1,r)   ;%floor(a1*sin(b1*(1:r)+c1) + a2*sin(b2*(1:r)+c2))+30;





















%% [Mehul] 13/March/2020 - OLD CODE:For testing purposes replaced with a new one
% % Img = imadjust(Img);
% [p, q, r] = size(Img);
% 
% %Img_vect = reshape(Img,m*n*l,1);
% %Img_vect = sort(Img_vect);
% 
% thresholds = zeros(r, 1);
% %size(Img_vect)
% %my_thresh = Img_vect(round(m*n*l*0.997),1)
% 
% for i = 1:r
%     % 'found local maxima within a box with size approximating the size of
%     % basal body 0.25x0.25x0.6 (um)^3'
%     % i.e. each plane taken to have thickness of 3~5 pixels
%     if i == 1
%         currPlanes = Img(:, :, i:i+2);
%         currBWPlanes = BW(:, :, i:i+2);
%         % row, col, depth coordinates of identified regional maxima
%         % note that 'find(X)' returns linear indices
%         [row, col, depth] =  ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%     elseif i == 2
%         currPlanes = Img(:, :, i-1:i+2);
%         currBWPlanes = BW(:, :, i-1:i+2);
%         [row, col, depth] = ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%     elseif i == r-1
%         currPlanes = Img(:, :, i-2:i+1);
%         currBWPlanes = BW(:, :, i-2:i+1);
%         [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%     elseif i==r
%         currPlanes=Img(:, :, i-2:i);
%         currBWPlanes=BW(:, :, i-2:i);
%         [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%     else
%         currPlanes=Img(:, :, i-2:i+2);
%         currBWPlanes=BW(:, :, i-2:i+2);
%         [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%     end
%     
%     currPlaneIntensity = zeros(length(row), 1);
%     for pt = 1:length(row)
%         currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
%     end
%     
%     %currMean = mean(currPlaneIntensity)
%     currMean = mean(currPlaneIntensity,'all');    
%     currStd = std(currPlaneIntensity);
%     currCv = currStd/currMean;
%     
%     % pg. 17 of Jingyi's paper: 'intensity threshold for given plane ...'
% %     currThreshold = (th - th*currCv)*currStd + currMean; % given in Chad's paper
% %     %((6?(Icv?6))×Is.d.)+Iavg
%     currThreshold = (th-(currCv*th))*currStd + currMean;
%     
% %     my_max = max(Img,[],'all');
%     mean(currPlanes);
%     thresholds(i) = currThreshold; %mean(mean(mean(currPlanes)))*50;
%     thresholds(i) = mean(currPlanes,'all')*5;
%     
%     
% end
% 
% % plot(thresholds)
% mtx = zeros(size(BW));
% mtx_temp = zeros(size(BW));
% % [cell1, ~] = identifyCell(Img);
% cell2 = calibrate_detect(Img);
% % imshow(cell)
% hold on
% sum(BW,'all')
% 
% 
% for i=1:p
%     for j=1:q
%         for k=1:r
%             % adaptive threshold used to filter local maxima
%             if BW(i, j, k) == 1 && Img(i, j, k) >= thresholds(k)
% %                 if cell1(i,j)==1 || cell2(i,j)==1
%                 mtx(i, j, k) = Img(i, j, k);
% %                 end
%             end
%         end
%     end
% end
% 
% 
% 
% [y, x, z] = ind2sub(size(mtx), find(mtx > 0));
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 







%     switch(i)
%         case 1
%             currPlanes = Img(:, :, i:i+2);
%             currBWPlanes = BW(:, :, i:i+2);
%             % row, col, depth coordinates of identified regional maxima
%             % note that 'find(X)' returns linear indices
%             [row, col, depth] =  ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%         case 2
%             currPlanes = Img(:, :, i-1:i+2);
%             currBWPlanes = BW(:, :, i-1:i+2);
%             [row, col, depth] = ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%         case r-1
%             currPlanes = Img(:, :, i-2:i+1);
%             currBWPlanes = BW(:, :, i-2:i+1);
%             [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%         case r
%             currPlanes=Img(:, :, i-2:i);
%             currBWPlanes=BW(:, :, i-2:i);
%             [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%         otherwise
%             currPlanes=Img(:, :, i-2:i+2);
%             currBWPlanes=BW(:, :, i-2:i+2);
%             [row, col, depth]= ind2sub(size(currBWPlanes), find(currBWPlanes > 0));
%     end
%     
%     currPlaneIntensity = zeros(length(row), 1);
%     for pt = 1:length(row)
%         currPlaneIntensity(pt) = currPlanes(row(pt), col(pt), depth(pt));
%     end
%     
%     %currMean = mean(currPlaneIntensity)
%     currMean = mean(mean(mean(currPlanes)));
%     mean(mean(mean(currBWPlanes)));
%     
%     currStd = std(currPlaneIntensity);
%     
%     currCv = currStd/currMean;
%     
%     % pg. 17 of Jingyi's paper: 'intensity threshold for given plane ...'
% %     currThreshold = (th - th*currCv)*currStd + currMean; % given in Chad's paper
%     %((6?(Icv?6))?Is.d.)+Iavg
%     currThreshold = (th-(currCv*th))*currStd + currMean;
%     
%     % thresholds(i) = currThreshold;
%     %prctile(Img,1,'all')
%     %thresholds(i) = prctile(Img,99,'all')
%     %my_thresh = graythresh(currPlanes)
%     %thresholds(i) = prctile(currPlaneIntensity,my_thresh,'all')
%     my_max = max(max(max(Img)));
%     %size(my_thresh)
%     %size(my_max)
%     %thresholds(i) = currMean + 2 * currStd
%     mean(currPlanes);
%     thresholds(i) = currThreshold;   %mean(mean(mean(currPlanes)))*50;
% %     thresholds(i) = mean(mean(mean(currPlanes)))*5;
%     %thresholds(i) = mean(mean(mean(currPlanes)))*10;
%     %thresholds(i) = max(max(max(currPlanes))) * mean(mean(mean(currPlanes)))
%     %thresholds(i,1) = currCv
%     %thresholds(i,2) = currMean;
%     %thresholds(i,3) = currStd
%     %thresholds(i,4) = currThreshold
%     
% end





% thresholds(i) = currThreshold;
    %prctile(Img,1,'all')
    %thresholds(i) = prctile(Img,99,'all')
    %my_thresh = graythresh(currPlanes)
    %thresholds(i) = prctile(currPlaneIntensity,my_thresh,'all')
    %size(my_thresh)
    %size(my_max)
    %thresholds(i) = currMean + 2 * currStd
    %thresholds(i) = mean(mean(mean(currPlanes)))*10;
    %thresholds(i) = max(max(max(currPlanes))) * mean(mean(mean(currPlanes)))
    %thresholds(i,1) = currCv
    %thresholds(i,2) = currMean;
    %thresholds(i,3) = currStd
    %thresholds(i,4) = currThreshold -  -