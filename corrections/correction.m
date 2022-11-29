function correction(I, rejection_threshold, visualization, minBBsInRow, minRowLength, resultPath, imageID, Iori)
    channel_count = 1;
    size(I);
    format shortg

    scale_xyz = [1, 1, 1];

    [cort_x, cort_y, cort_z, antPole, postPole, order2OA, pot_oa, BW, I2] = getBBIdx_4correction(I, rejection_threshold, scale_xyz, resultPath, imageID, Iori);
    % [cort_x, cort_y, cort_z, oa_x, oa_y, oa_z] = getBBIdx2(I, [0.108, 0.108, 0.1]);
    oa_x = pot_oa(:, 1);
    oa_y = pot_oa(:, 2);
    oa_z = pot_oa(:, 3);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % repeated calculation of antPole and postPole within ... remove
    % redundancy later
    % note that we are going to calculate final_traceback based on
    % cortical BBs and are excluded OA BBs for now.
    for c = 1: 20
        cort_x = cort_x*.75;
        cort_y = cort_y*.75;
        [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows4(cort_x, cort_y, cort_z, antPole, postPole, minBBsInRow, minRowLength, scale_xyz);
        cort_x = cort_x/.75;
        cort_y = cort_y/.75;
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

    BBstat2(I2, BW,cort_x, cort_y, cort_z, resultPath, Iori)
end