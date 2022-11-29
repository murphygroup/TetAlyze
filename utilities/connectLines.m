% -------------------------------------------------------------------------
% [Paul] 9/5/20
% -------------------------------------------------------------------------


function [label, updatedTraceback] = connectLines(mtx, label, dist2Ant, dist2Post, dist2All, x, y, z, antPole, postPole)
% k: row index, i bb index
numBBs_row = sum(mtx ~= 0, 2);
tails = zeros(length(mtx),1);
[d1, ~] = size(mtx);
for i = 1:d1
    tails(i) = mtx(i, numBBs_row(i));
end
k = 1;
while true
    heads = mtx(:, 1);
    numBBs_row = sum(mtx ~= 0, 2);
    [d1, ~] = size(mtx);
    tails = zeros(d1,1);
    for j = 1:d1
        num = numBBs_row(j);
        tails(j) = mtx(j, num);
    end
    sig = false;
    head = heads(k);
    dists = dist2All(head, :);
    [~, I] = sort(dists, 'ascend'); 
    for i = 2:9
        if ~ismember(I(i),tails)
            continue;
        else
            rowNum1 = numBBs_row(k);
            k2 = find(tails == I(i));
            rowNum2 = numBBs_row(k2);
            currRow1 = mtx(k, 1:rowNum1);
            currRow2 = mtx(k2, 1:rowNum2);
            pt1 = [mean(x(currRow1)), mean(y(currRow1)), mean(z(currRow1))];
            pt2 = [mean(x(currRow2)), mean(y(currRow2)), mean(z(currRow2))];
            vec1 = pt1 - pt2;
            vec2 = antPole-postPole;
            angle = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
            angle = angle/pi*180;
            
            bb = [x(mtx(k, 2)), y(mtx(k, 2)), z(mtx(k, 2))];
            bb2 = [x(mtx(k2, rowNum2-1)), y(mtx(k2, rowNum2-1)), z(mtx(k2, rowNum2-1))];
            bbHead3d = [x(mtx(k, 1)), y(mtx(k, 1)), z(mtx(k, 1))];
            bbTail3d = [x(mtx(k2, rowNum2)), y(mtx(k2, rowNum2)), z(mtx(k2, rowNum2))];
            vec1 = bbHead3d - bb;
            vec2 = bbTail3d - bbHead3d;
            vec3 = bbTail3d - bb2;
            angle1 = acos(dot(vec1, vec2)/(norm(vec1)*norm(vec2)));
            angle1 = angle1/pi*180;
            angle2 = acos(dot(vec3, vec2)/(norm(vec3)*norm(vec2)));
            angle2 = angle2/pi*180;
            angle3 = fourptsangle(bbHead3d, bbTail3d, antPole, postPole);
            angle3 = angle3/pi*180;
            
            row = zeros(rowNum1+rowNum2, 1);
            row(1:rowNum1) = mtx(k, 1:rowNum1);
            row(rowNum1+1:end) = mtx(k2, 1:rowNum2);
            row = sortByRow(row', 1, dist2Ant);
            if angle2 > 120 && angle1 < 60 && angle3 < 10  && checkOneRow(row, x, y, z)
                % merge row k, k2
                mtx(k, rowNum1+1:rowNum1+rowNum2) = mtx(k2,1:rowNum2);
                % update 'label' vector
                label(label == k2) = k;
                mtx = sortByRow(mtx, k, dist2Ant);
                mtx(k2,:)=[];
                sig = true;
                break;
            end
        end
    end
    if ~sig
        k = k + 1;
    end
    [d1, ~] = size(mtx);
    if k > d1
        break
    end
end
    
updatedTraceback = mtx;
end