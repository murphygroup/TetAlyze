% New version
function [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows2(x, y, z, oa_x, oa_y, oa_z, antPole, postPole,~,scale_xyz)

num_BBs = length(x);
all_pts = [x,y,z].*scale_xyz;
dist2Ant = vecnorm(all_pts - antPole.*scale_xyz,2,2);
dist2Post = vecnorm(all_pts - postPole.*scale_xyz,2,2);

dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.7,0.3,0, 0.0]);
[startPt, endPt] = findLink(x, y, z, dist2All, antPole, postPole, dist2Ant, dist2Post);
% matrix showing ciliary rows
OriginalBBsMatrix = link(startPt, endPt);
[d1, ~] = size(OriginalBBsMatrix);
for k=1:d1
    OriginalBBsMatrix = sortByRow(OriginalBBsMatrix, k, dist2Ant);
end
updatedTraceback = OriginalBBsMatrix;
NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = find(NPR <= 3);
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

% x = x/.5;
% y = y/.5;
% figure(2)
% visualize_bbs(x, y, z, 3*[], 3*[],[], updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);

changed = true;
ant2post = postPole - antPole;
counter = 0;
while changed && counter < 100000
    changed = false;
    counter = counter + 1;
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);

    [d1, ~] = size(updatedTraceback);
    NPR = sum(updatedTraceback ~= 0, 2);
    heads = updatedTraceback(:, 1);
    tails = zeros(d1, 1);
    for i = 1:d1
        tails(i) = updatedTraceback(i, NPR(i));
    end

    L = 2 * length(tails) + length(withoutLabel);
    candidates = zeros(L, 1);
    candidates(1:length(tails)) = tails;
    candidates(length(tails)+1:2*length(tails)) = heads;
    candidates(2*length(tails)+1:length(candidates)) = withoutLabel;
    distances = dist2All(candidates, candidates);
    [startPt, endPt] = findClosest(candidates, distances, heads, tails);
    for i = 1:length(startPt)
        if (ismember(startPt(i), heads) && ismember(endPt(i), tails))
            p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
            if angleWithPoles < 50/180*pi || (angleWithPoles < 80/180*pi && (dist2Ant(startPt(i)) < 10 || (dist2Post(startPt(i)) < 10)))
                k1 = find(heads == startPt(i));
                k2 = find(tails == endPt(i));
                if k1 ~= k2
                    updatedTraceback = connect1Line(updatedTraceback, k1, k2, NPR(k1), NPR(k2), dist2Ant, dist2All);
                    changed = true;
                    break
                end
            end
        elseif (ismember(startPt(i), tails) && ismember(endPt(i), heads))
            p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
            if angleWithPoles < 50/180*pi || (angleWithPoles < 80/180*pi && (dist2Ant(startPt(i)) < 10 || (dist2Post(startPt(i)) < 10)))
                k1 = find(tails == startPt(i));
                k2 = find(heads == endPt(i));
                if k1 ~= k2
                    updatedTraceback = connect1Line(updatedTraceback, k2, k1, NPR(k2), NPR(k1), dist2Ant, dist2All);
                    changed = true;
                    break
                end
            end

        elseif (ismember(startPt(i), heads) && ismember(endPt(i), withoutLabel))
            k = find(heads == startPt(i));
            vec1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            vec2 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - ...
                   [x(updatedTraceback(k, 2)), y(updatedTraceback(k, 2)), z(updatedTraceback(k, 2))];
            angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
            angleLocal = angleLocal*180/pi;
            if angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
                updatedTraceback(k, NPR(k)+1) = endPt(i);
                updatedTraceback = sortByRow(updatedTraceback, k, dist2Ant);
                changed = true;
                break
            end
        elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), heads))
            k = find(heads == endPt(i));
            vec1 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - [x(startPt(i)), y(startPt(i)), z(startPt(i))];
            vec2 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - ...
                   [x(updatedTraceback(k, 2)), y(updatedTraceback(k, 2)), z(updatedTraceback(k, 2))];
            angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
            angleLocal = angleLocal*180/pi;
            if angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
                updatedTraceback(k, NPR(k)+1) = startPt(i);
                updatedTraceback = sortByRow(updatedTraceback, k, dist2Ant);
                changed = true;
                break
            end
        elseif (ismember(startPt(i), tails) && ismember(endPt(i), withoutLabel))
            k2 = find(tails == startPt(i));
            vec1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            vec2 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - ...
                   [x(updatedTraceback(k2, NPR(k2)-1)), y(updatedTraceback(k2, NPR(k2)-1)), z(updatedTraceback(k2, NPR(k2)-1))];
            angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
            angleLocal = angleLocal*180/pi;
            if angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
                updatedTraceback(k2, NPR(k2)+1) = endPt(i);
                updatedTraceback = sortByRow(updatedTraceback, k2, dist2Ant);
                changed = true;
                break
            end
        elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), tails))
            k2 = find(tails == endPt(i));
            vec1 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - [x(startPt(i)), y(startPt(i)), z(startPt(i))];
            vec2 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - ...
                   [x(updatedTraceback(k2, NPR(k2)-1)), y(updatedTraceback(k2, NPR(k2)-1)), z(updatedTraceback(k2, NPR(k2)-1))];
            angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
            angleLocal = angleLocal*180/pi;
            if angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
                updatedTraceback(k2, NPR(k2)+1) = startPt(i);
                updatedTraceback = sortByRow(updatedTraceback, k2, dist2Ant);
                changed = true;
                break
            end
        elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), withoutLabel))
            p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
            if angleWithPoles < 35/180*pi || (angleWithPoles < 80/180*pi && (dist2Ant(startPt(i)) < 10 || (dist2Post(startPt(i)) < 10)))
                updatedTraceback(end+1,1:2) = [startPt(i), endPt(i)];
                [d1, ~] = size(updatedTraceback);
                updatedTraceback = sortByRow(updatedTraceback, d1, dist2Ant);
                changed = true;
                break
            end
        end

    end
end
NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = find(NPR <= 4);
updatedTraceback(deleteRows, :)=[];


% updatedTraceback = checkRows(updatedTraceback, x, y, z);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

% x = x/.5;
% y = y/.5;
% figure(2)
% visualize_bbs(x, y, z, 3*[], 3*[],[], updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
minBBsInRow = 2;
for i=1:10
    [~, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, ...
        antPole, postPole, minBBsInRow, scale_xyz, dist2Post);
    withLabel = intersect(1:num_BBs, updatedTraceback);

    [~, updatedTraceback] = connectLines(updatedTraceback, withLabel, dist2Ant, dist2Post, dist2All, x, y, z, antPole, postPole);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
end
% x = x/.5;
% y = y/.5;
% figure(2)
% visualize_bbs(x, y, z, 3*[], 3*[],[], updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = find(NPR <= 3);
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

end