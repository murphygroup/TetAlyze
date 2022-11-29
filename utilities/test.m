clear all
tempPath = 'D:\Tet generative model\generated data';
for q=0:180
    close all hidden
    filename = sprintf('WTv2 %d.txt', q);
    cort_x = table2array(readtable(fullfile(tempPath, sprintf('xWTv2_%d.csv', q))));
    cort_y = table2array(readtable(fullfile(tempPath, sprintf('yWTv2_%d.csv', q))));
    cort_z = table2array(readtable(fullfile(tempPath, sprintf('zWTv2_%d.csv', q))));
    resultFolderPath = 'D:\Tet generative model\test';
    updatedTraceback = table2array(readtable(fullfile(tempPath, sprintf('WT_alignment_%d.csv', q))));

    rejection_threshold = 3;  % 0.5 (strong rejection) - 3 (weak rejection) for WT or big1 without concaves, 3-7 for Big 1 with concaves
    minBBsInRow = 3; % Rows with fewer BBs than this value will be discarded.
    minRowLength = 0; % Threshold for BB row length.
    imageID = '/';
%     antPole = [0 0 50];
%     postPole = [0 0 0];
%     scale_xyz = [1, 1, 1];
    for c = 1: 1
        tails = [];
%         cort_x = cort_x*.75;
%         cort_y = cort_y*.75;
%         [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows4(cort_x, cort_y, cort_z, antPole, postPole, minBBsInRow, minRowLength, scale_xyz);
%         cort_x = cort_x/.75;
%         cort_y = cort_y/.75;
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
        [cort_x, cort_y, cort_z, oa_x, oa_y, oa_z] = newCoor(cort_x, cort_y, cort_z, [0], [0], [0], antPole, postPole);
        antPole = [0 0 norm(antPole - postPole)];
        postPole = [0 0 0];
    end

    antIdx = [];
    postIdx = [];
    midIdx = [];
    a = max(cort_z) ;
    b = min(cort_z);
    for i=1:length(cort_x)
        if cort_z(i) >= a - 0.25*(a-b)
            antIdx(end+1) = i;
        elseif cort_z(i) <= b + 0.25*(a-b)
            postIdx(end+1) = i;
        else
            midIdx(end+1) = i;
        end
    end
    surf(cort_x, cort_x, cort_y, cort_z, antIdx, postIdx, midIdx, 'D:\Tet generative model\test', filename)

    [num_rows, ~] = size(updatedTraceback);

    xy_centroids = zeros(num_rows, 1);
    for i = 1:num_rows
        row = updatedTraceback(i, :);
        row = row(row ~= 0);
        xy_centroids(i) = atan2(mean(cort_y(row)), mean(cort_x(row)));
    end
    [~, order] = sort(xy_centroids, 'ascend');

    % for i=1:length(order2OA)
    %     if ismember(order2OA(i), withLabel)
    %         [startRow,~] = find(updatedTraceback==order2OA(i));
    %         break
    %     end
    % end
    startRow = 1;
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

%     figure(2)
%     hFig = figure(2);
    % set(gcf,'position',get(0,'ScreenSize'))
%     visualize_bbs(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z,  updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
    % saveas(hFig, join([resultPath, imageID, 'Alignment.m']));
%     prompt = 'Do you want to correct the row index? If yes, enter the current index of the real No.1 row and press Enter. Or, directly press Enter.\n';
%     str = input(prompt,'s');
%     if ~isempty(str)
%     idx = str2num(str);
    newOrder2 = zeros(length(order), 1);
    for k=1:length(newOrder)
        id = idx+k-1;
        if id > length(newOrder)
            id = id - length(newOrder);
        end
        newOrder2(k) = id;

    end
    updatedTraceback = updatedTraceback(newOrder2, :);
%     end

    vec = analysis(updatedTraceback, cort_x, cort_y, cort_z, antPole, postPole, 'D:\Tet generative model\test', imageID, filename);
    f2 = join([resultFolderPath, '/data4PCA.csv']);
    writematrix(vec,f2);
%     if size(findobj(hFig))>0
%         close(hFig)
%     end
end