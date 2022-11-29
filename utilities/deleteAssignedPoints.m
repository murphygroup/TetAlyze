% -------------------------------------------------------------------------
% [Paul] 9/15/20
% -------------------------------------------------------------------------


function mtx = deleteAssignedPoints(num, mtx, scale, x, y, z, antPole, postPole)
    numBBs_row = sum(mtx ~= 0, 2);
    [d1, ~] = size(mtx);
    selectBBs = zeros(1,d1);
    for i = 1:d1
        row = mtx(i,1:numBBs_row(i));
        distance = zeros(1,length(row));
        ptMean = [mean(x(row)), mean(y(row)), mean(z(row))];
        for j =1:length(row)
            idx = row(j);
            pt = [x(idx), y(idx), z(idx)];
            distance(j) = distance_pt2plane(pt, scale, ptMean, antPole, postPole);
        end
        [~,idx] = max(distance);
%         idx = randi(numBBs_row(i), 1);
        selectBBs(i) = mtx(i, idx);
    end

    for i = 1:length(selectBBs)
    %     r = rand;
        if  numBBs_row(i)>3  % delete point
            BB = selectBBs(i);
            idx = find(mtx(i, 1:numBBs_row(i))==BB);
            if idx == 1
                mtx(i,1:numBBs_row(i)-1) = mtx(i, 2:numBBs_row(i));
                mtx(i,numBBs_row(i)) = 0;
            elseif idx == numBBs_row(i)
                mtx(i,numBBs_row(i)) = 0;
            elseif idx == 2
                mtx(i,1:numBBs_row(i)-2) = mtx(i, 3:numBBs_row(i));
                mtx(i,numBBs_row(i)-1) = 0;
                mtx(i,numBBs_row(i)) = 0;
            elseif idx == numBBs_row(i)-1
                mtx(i,numBBs_row(i)-1:numBBs_row(i)) = 0;
            else
                mtx(end+1,1:numBBs_row(i)-idx) = mtx(i, idx+1:numBBs_row(i));
                mtx(i,idx+1:numBBs_row(i)) = 0;
            end
        end
    end
end


        