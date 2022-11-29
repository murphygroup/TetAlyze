% -------------------------------------------------------------------------
% [Ben] 05/28/18
% -------------------------------------------------------------------------

function vec = train_model(I, thresh_ratio, rejection_threshold, visualization, minBBsInRow, minRowLength, resultPath, imageID, Iori)
channel_count = 1;
size(I);
format short

scale_xyz = [1, 1, 1];

[cort_x, cort_y, cort_z, antPole, postPole, order2OA, pot_oa, BW, I2, Iori] = getBBIdx4(I, thresh_ratio, rejection_threshold, scale_xyz, resultPath, imageID, Iori);
% [cort_x, cort_y, cort_z, oa_x, oa_y, oa_z] = getBBIdx2(I, [0.108, 0.108, 0.1]);




oa_x = pot_oa(:, 1);
oa_y = pot_oa(:, 2);
oa_z = pot_oa(:, 3);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% repeated calculation of antPole and postPole within ... remove
% redundancy later
% note that we are going to calculate final_traceback based on
% cortical BBs and are excluded OA BBs for now.
for c = 1: 30
%     cort_x = cort_x*.75;
%     cort_y = cort_y*.75;
    [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows3(cort_x, cort_y, cort_z, antPole, postPole, minBBsInRow, minRowLength, scale_xyz);
%     cort_x = cort_x/.75;
%     cort_y = cort_y/.75;
    [d1, ~] = size(updatedTraceback);
    NPR = sum(updatedTraceback ~= 0, 2);
    heads = updatedTraceback(:, 1);
    for i = 1:d1
        tails(i) = updatedTraceback(i, NPR(i));
    end
    idx1 = [];
    idx2 = [];
    headBB = [cort_x(heads), cort_y(heads), cort_z(heads)];
    meanH = mean(headBB(:, 3));
    for i = 1:d1
        if headBB(i, 3) < meanH
            idx1(end+1) = i;
        else
            idx2(end+1) = i;
        end
    end
    if length(idx1) > length(idx2)
        headBB = headBB(idx1, :);
    else
        headBB = headBB(idx2, :);
    end
    meanH = mean(headBB, 1);
    idx1 = [];
    idx2 = [];
    tailBB = [cort_x(tails), cort_y(tails), cort_z(tails)];
    meanT = mean(tailBB(:, 3));
    for i = 1:d1
        if tailBB(i, 3) < meanT
            idx1(end+1) = i;
        else
            idx2(end+1) = i;
        end
    end
    if length(idx1) > length(idx2)
        tailBB = tailBB(idx1, :);
    else
        tailBB = tailBB(idx2, :);
    end
    meanT = mean(tailBB, 1);
    meanT(3) = meanT(3) + 2;
    meanH(3) = meanH(3) -2;
    % meanDist = stats(updatedTraceback, scoreMtx, cort_x, cort_y, cort_z, antPole, postPole, scale_xyz);
    antPole = meanT;
    postPole = meanH;
    [cort_x, cort_y, cort_z, oa_x, oa_y, oa_z] = newCoor(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, antPole, postPole);
    antPole = [0 0 norm(antPole - postPole)];
    postPole = [0 0 0];
    
end

BBstat(I2, BW,cort_x, cort_y, cort_z, resultPath, Iori)
[num_rows, ~] = size(updatedTraceback);

xy_centroids = zeros(num_rows, 1);
for i = 1:num_rows
    row = updatedTraceback(i, :);
    row = row(row ~= 0);
    xy_centroids(i) = atan2(mean(cort_y(row)), mean(cort_x(row)));
end
[~, order] = sort(xy_centroids, 'ascend');

for i=1:length(order2OA)
    if ismember(order2OA(i), withLabel)
        [startRow,~] = find(updatedTraceback==order2OA(i));
        break
    end
end
newOrder = zeros(length(order), 1);

idx = find(order == startRow);
for k=1:length(order)
    newOrder(k) = order(idx);
    idx = idx + 1;
    if idx > length(order)
        idx = idx - length(order);
    end
end
updatedTraceback = updatedTraceback(newOrder, :);
% 
if visualization
    hFig = figure(2);
    % set(gcf,'position',get(0,'ScreenSize'))
    visualize_bbs(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z,  updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
    % saveas(hFig, join([resultPath, imageID, 'Alignment.m']));
    prompt = 'Do you want to correct the row index? If yes, enter the current index of the real No.1 row and press Enter. Or, directly press Enter.\n';
    str = input(prompt,'s');
    if ~isempty(str)
        idx = str2num(str);
        if idx < 0
            temp = antPole;
            antPole = postPole;
            postPole = temp;
            a = length(oa_x);
            cort_x(end+1:end+a) = pot_oa(:, 1);
            cort_y(end+1:end+a) = pot_oa(:, 2);
            cort_z(end+1:end+a) = pot_oa(:, 3);

            [cort_x, cort_y, cort_z] = newCoorWithoutOA(cort_x, cort_y, cort_z, antPole, postPole);

            pot_oa(:, 1) = cort_x(end-a+1:end);
            pot_oa(:, 2) = cort_y(end-a+1:end);
            pot_oa(:, 3) = cort_z(end-a+1:end);
            cort_x(end-a+1:end) = [];
            cort_y(end-a+1:end) = [];
            cort_z(end-a+1:end) = [];

            [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows3(cort_x, cort_y, cort_z, antPole, postPole, minBBsInRow, minRowLength, scale_xyz);

            [d1, ~] = size(updatedTraceback);
            NPR = sum(updatedTraceback ~= 0, 2);
            close(hFig)


            figure(2)
            hFig = figure(2);
            % set(gcf,'position',get(0,'ScreenSize'))
            visualize_bbs(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z,  updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
            prompt = 'Do you want to correct the row index? If yes, enter the current index of the real No.1 row and press Enter. Or, directly press Enter.\n';
            str = input(prompt,'s');
            idx = str2num(str);
        end
        newOrder2 = zeros(length(order), 1);
        for k=1:length(newOrder)
            id = idx+k-1;
            if id > length(newOrder)
                id = id - length(newOrder);
            end
            newOrder2(k) = id;

        end
        updatedTraceback = updatedTraceback(newOrder2, :);
    end
    writematrix(updatedTraceback, join([resultPath, 'Alignment.csv']));

    if size(findobj(hFig))>0
        close(hFig)
    end
end

if visualization
    figure(2)
    hFig = figure(2);
%     set(gcf,'position',get(0,'ScreenSize'))
    visualize_bbs(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
    saveas(hFig, join([resultPath, 'Alignment.m']));
end

immBBnum = intensityAnalysis(cort_x, cort_y, cort_z, BW, I2, updatedTraceback(1,:), resultPath);
vec = analysis(updatedTraceback, cort_x, cort_y, cort_z, antPole, postPole, resultPath, imageID, 'Summary.txt', immBBnum);

end