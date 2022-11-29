% -------------------------------------------------------------------------
% [Ben] 2/1/18
% Returns an updated traceback matrix, together with the corresponding
% label vector. The traceback matrix consists of sequences of BBs (padded
% with trailing 0's) in each row, tracing out a ciliary row. The label
% vector indicates which BBs belong to the same ciliary row. The updated
% traceback matrix may have fewer rows than the original traceback matrix.
% In effect, this function deletes ciliary rows with too few BBs and
% reassigns the BBs from deleted ciliary rows to other ciliary rows.
% -------------------------------------------------------------------------


function [iniLabel, updatedTraceback] = getIniLabel(mtx, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, ...
    antPole, postPole, minBBsInRow, scale_xyz, dist2Post)
% we assume the following invariants on the input arguments:
% 1) the BB indices in each row of mtx are already sorted in order of
% decreasing distance from the anterior pole.

[numRows, ~] = size(mtx);
CRs = 1:numRows;
row_counts = sum(mtx ~= 0, 2);
legalRows = CRs(row_counts >= 2);
legal_row_counts = row_counts(legalRows);
numLegalRows = length(legalRows);

% numRows is the no. of ciliary rows detected by findLink.m and link.m
numBBs = length(x);
label = zeros(numBBs, 1);

for i = 1:numLegalRows
    row = legalRows(i);
    bbs = mtx(row, 1:row_counts(row));
    label(bbs) = i;
end



% -------------------------------------------------------------------------
% [Chongming] 8/20/20
% consider the dihedral angle of plane (pt1, antPole, postPole) and (pt2, antPole, postPole)
% -------------------------------------------------------------------------
% k: row index, i bb index
k = 1;
while k <= numRows
    sig = false;
    numBBs_row = sum(mtx(k, :) ~= 0, 2);
    bbHead = mtx(k, 1);
    added = 0;
% -------------------------------------------------------------------------
% search from head
    dists = dist2All(bbHead, :);
    [~, I] = sort(dists, 'ascend'); 
    candidateIdx = [];
    for i = 2:16
        if ismember(I(i),withLabel)
            continue;
        else
            candidateIdx(end+1) = I(i);
        end
        if length(candidateIdx) > 2
            break
        end
    end
    minIdx = -1;
    minAngle = 400;
    for i=1:length(candidateIdx)
        idx = candidateIdx(i);
        bb = [x(idx), y(idx), z(idx)];
        bbHead3d = [x(bbHead), y(bbHead), z(bbHead)];  
        bb2 = [x(mtx(k, 2)), y(mtx(k, 2)), z(mtx(k, 2))];
        localVec1 = bbHead3d - bb2;
        localVec2 = bbHead3d - bb;
        ant2post = antPole - postPole;
        p22p1 = bbHead3d - bb;
        angleLocal = acos(dot(localVec1, localVec2)/(norm(localVec1)*norm(localVec2)));
        angleLocal = angleLocal/pi*180;
        angleWithPoles =  fourptsangle(bb, bbHead3d, antPole, postPole);
%         acos(abs(dot(ant2post, p22p1))/(norm(ant2post)*norm(p22p1)));
        angleWithPoles = angleWithPoles/pi*180;
        angleLocal3 = 0;
        if numBBs_row > 2
            bb3 = [x(mtx(k, 3)), y(mtx(k, 3)), z(mtx(k, 3))];
            localVec1 = bb2 - bb3;
            localVec2 = bb2 - bbHead3d;
            angleLocal2 = acos(dot(localVec1, localVec2)/(norm(localVec1)*norm(localVec2)));
            angleLocal2 = angleLocal2/pi*180;
            angleLocal3 = abs(angleLocal-angleLocal2);
        end

        if angleLocal3 < 45 && angleLocal > 120 && angleWithPoles < 45
            if true
                legal_row_counts(k) = legal_row_counts(k) + 1;
                mtx(k, legal_row_counts(k)) = idx;
                % update 'label' vector
                label(idx) = k;
                mtx = sortByRow(mtx, k, dist2Ant);
                withLabel(end+1) = idx;
                sig = true;
                added = added + 1;
            end
        end
        if sig
            numBBs_row = sum(mtx(k, :) ~= 0, 2);
            bbHead = mtx(k, 1);
        end
        if added > 2
            break
        end
    end
    if ~sig || added > 2
        k = k + 1;
    end
    [numRows, ~] = size(mtx);
end
k = 1;
while k <= numRows
    sig = false;
% -------------------------------------------------------------------------
% search from tail
    numBBs_row = sum(mtx(k, :) ~= 0, 2);
    bbTail = mtx(k, numBBs_row);
    dists = dist2All(bbTail, :);
    [~, I] = sort(dists, 'ascend'); 

    candidateIdx = [];
    for i = 2:16
        if ismember(I(i),withLabel)
            continue;
        else
            candidateIdx(end+1) = I(i);
        end
        if length(candidateIdx) > 2
            break
        end
    end
    minIdx = -1;
    minAngle = 400;
    for i=1:length(candidateIdx)
        idx = candidateIdx(i);
        bb = [x(idx), y(idx), z(idx)];
        bbTail3d = [x(bbTail), y(bbTail), z(bbTail)];
        bb2 = [x(mtx(k, numBBs_row-1)), y(mtx(k, numBBs_row-1)), z(mtx(k, numBBs_row-1))];
        localVec1 = bbTail3d - bb2;
        localVec2 = bbTail3d - bb;
        ant2post = antPole - postPole;
        p22p1 = bbTail3d - bb;
        angleLocal = acos(dot(localVec1, localVec2)/(norm(localVec1)*norm(localVec2)));
        angleLocal = angleLocal/pi*180;
        angleWithPoles = fourptsangle(bb, bbTail3d, antPole, postPole);
%         acos(abs(dot(ant2post, p22p1))/(norm(ant2post)*norm(p22p1)));
        angleWithPoles = angleWithPoles/pi*180;
        angle3 = fourptsangle(bb, bbTail3d, antPole, postPole);
        angle3 = angle3/pi*180;
        angleLocal3 = 0;
        if numBBs_row > 2
            bb3 = [x(mtx(k, numBBs_row-2)), y(mtx(k, numBBs_row-2)), z(mtx(k, numBBs_row-2))];
            localVec1 = bb2 - bb3;
            localVec2 = bb2 - bbTail3d;
            angleLocal2 = acos(dot(localVec1, localVec2)/(norm(localVec1)*norm(localVec2)));
            angleLocal2 = angleLocal2/pi*180;
            angleLocal3 = abs(angleLocal-angleLocal2);
        end

        if  angleLocal3 < 45 && angleLocal > 120 && angleWithPoles < 45 
            if true
                legal_row_counts(k) = legal_row_counts(k) + 1;
                mtx(k, legal_row_counts(k)) = idx;
                % update 'label' vector
                label(idx) = k;
                mtx = sortByRow(mtx, k, dist2Ant);
                withLabel(end+1) = idx;
                sig = true;
                added = added + 1;
            end
        end
        if sig
            numBBs_row = sum(mtx(k, :) ~= 0, 2);
            bbTail = mtx(k, numBBs_row);
        end
        if added > 2
            break
        end
    end
    if ~sig || added > 2
        k = k + 1;
    end
    [numRows, ~] = size(mtx);
end
    
   
updatedTraceback = mtx;

    
iniLabel = label;
end