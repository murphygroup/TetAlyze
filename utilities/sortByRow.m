function oriMtx = sortByRow(oriMtx, row, dist2Ant)
% we assume the following invariants on the input arguments:
% 1) oriMtx does not have any trailing columns of all zeros on the right

row_counts = sum(oriMtx ~= 0, 2);
numBBs_row = row_counts(row);
bbs = oriMtx(row, 1:numBBs_row);
[~, I] = sort(dist2Ant(bbs), 'descend');
oriMtx(row, 1:numBBs_row) = bbs(I);
end