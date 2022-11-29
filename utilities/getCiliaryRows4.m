% New version
function [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows4(x, y, z, antPole, postPole, minBBsInRow, minRowLength, scale_xyz, minBBsInRow2)

num_BBs = length(x);
all_pts = [x,y,z].*scale_xyz;
dist2Ant = vecnorm(all_pts - antPole.*scale_xyz,2,2);
dist2Post = vecnorm(all_pts - postPole.*scale_xyz,2,2);

% dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.8, 0.2, 0, 0.0]);
dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.6, 0.4, 0, 0.0]);
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


% changed = true;
ant2post = postPole - antPole;
% counter = 0;
% while changed && counter < 10000
%     changed = false;
%     counter = counter + 1;
%     withLabel = intersect(1:num_BBs, updatedTraceback);
%     withoutLabel = setdiff(1:num_BBs, withLabel);
% 
%     [d1, ~] = size(updatedTraceback);
%     NPR = sum(updatedTraceback ~= 0, 2);
%     heads = updatedTraceback(:, 1);
%     tails = zeros(d1, 1);
%     for i = 1:d1
%         tails(i) = updatedTraceback(i, NPR(i));
%     end
% 
%     L = 2 * length(tails) + length(withoutLabel);
%     candidates = zeros(L, 1);
%     candidates(1:length(tails)) = tails;
%     candidates(length(tails)+1:2*length(tails)) = heads;
%     candidates(2*length(tails)+1:length(candidates)) = withoutLabel;
%     distances = dist2All(candidates, candidates);
%     [startPt, endPt] = findClosest(candidates, distances, heads, tails);
%     for i = 1:length(startPt)
%         if (ismember(startPt(i), heads) && ismember(endPt(i), tails))
%             p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
%             angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
%             if true % angleWithPoles < 50/180*pi || (angleWithPoles < 80/180*pi && (dist2Ant(startPt(i)) < 10 || (dist2Post(startPt(i)) < 10)))
%                 k1 = find(heads == startPt(i));
%                 k2 = find(tails == endPt(i));
%                 row1 = updatedTraceback(k1, 1:NPR(k1));
%                 row2 = updatedTraceback(k2, 1:NPR(k2));
%                 zRange1 = min(max(z(row1)) - min(z(row1)), max(z(row2)) - min(z(row2)));
%                 zRange2 = max(min(max(z(row1)), max(z(row2))) - max(min(z(row1)), min(z(row2))), 0);
%                 row = zeros(NPR(k1)+NPR(k2), 1);
%                 row(1:NPR(k1)) = row1;
%                 row(NPR(k1)+1:end) = row2;
%                 row = sortByRow(row', 1, dist2Ant);
%                 if k1 ~= k2 && checkOneRow(row, x, y, z)%&& (zRange2/(zRange1+0.00001) < 0.5)
%                     updatedTraceback = connect1Line(updatedTraceback, k1, k2, NPR(k1), NPR(k2), dist2Ant, dist2All);
%                     changed = true;
%                     break
%                 end
%             end
% 
% 
%         elseif (ismember(startPt(i), heads) && ismember(endPt(i), withoutLabel))
%             k = find(heads == startPt(i));
%             vec1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
%             vec2 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - ...
%                    [x(updatedTraceback(k, 2)), y(updatedTraceback(k, 2)), z(updatedTraceback(k, 2))];
%             angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
%             angleLocal = angleLocal*180/pi;
%             if true %angleLocal < 50 
%                 updatedTraceback(k, NPR(k)+1) = endPt(i);
%                 updatedTraceback = sortByRow(updatedTraceback, k, dist2Ant);
%                 changed = true;
%                 break
%             end
%         elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), heads))
%             k = find(heads == endPt(i));
%             vec1 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - [x(startPt(i)), y(startPt(i)), z(startPt(i))];
%             vec2 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - ...
%                    [x(updatedTraceback(k, 2)), y(updatedTraceback(k, 2)), z(updatedTraceback(k, 2))];
%             angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
%             angleLocal = angleLocal*180/pi;
%             if true %angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
%                 updatedTraceback(k, NPR(k)+1) = startPt(i);
%                 updatedTraceback = sortByRow(updatedTraceback, k, dist2Ant);
%                 changed = true;
%                 break
%             end
%         elseif (ismember(startPt(i), tails) && ismember(endPt(i), withoutLabel))
%             k2 = find(tails == startPt(i));
%             vec1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
%             vec2 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - ...
%                    [x(updatedTraceback(k2, NPR(k2)-1)), y(updatedTraceback(k2, NPR(k2)-1)), z(updatedTraceback(k2, NPR(k2)-1))];
%             angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
%             angleLocal = angleLocal*180/pi;
%             if true %angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
%                 updatedTraceback(k2, NPR(k2)+1) = endPt(i);
%                 updatedTraceback = sortByRow(updatedTraceback, k2, dist2Ant);
%                 changed = true;
%                 break
%             end
%         elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), tails))
%             k2 = find(tails == endPt(i));
%             vec1 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - [x(startPt(i)), y(startPt(i)), z(startPt(i))];
%             vec2 = [x(endPt(i)), y(endPt(i)), z(endPt(i))] - ...
%                    [x(updatedTraceback(k2, NPR(k2)-1)), y(updatedTraceback(k2, NPR(k2)-1)), z(updatedTraceback(k2, NPR(k2)-1))];
%             angleLocal = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
%             angleLocal = angleLocal*180/pi;
%             if true %angleLocal < 50 || (angleLocal < 80 && (dist2Ant(endPt(i)) < 10 || (dist2Post(endPt(i)) < 10)))
%                 updatedTraceback(k2, NPR(k2)+1) = startPt(i);
%                 updatedTraceback = sortByRow(updatedTraceback, k2, dist2Ant);
%                 changed = true;
%                 break
%             end
%         elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), withoutLabel))
%             p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
%             angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
%             if true %angleWithPoles < 70/180*pi || (angleWithPoles < 90/180*pi && (dist2Ant(startPt(i)) < 10 || (dist2Post(startPt(i)) < 10)))
%                 updatedTraceback(end+1,1:2) = [startPt(i), endPt(i)];
%                 [d1, ~] = size(updatedTraceback);
%                 updatedTraceback = sortByRow(updatedTraceback, d1, dist2Ant);
%                 changed = true;
%                 break
%             end
%         elseif (ismember(startPt(i), withoutLabel) && ismember(endPt(i), withoutLabel))
%             p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
% %             angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
%             if vecnorm(p22p1, 2) < 3
%                 updatedTraceback(end+1,1:2) = [startPt(i), endPt(i)];
%                 [d1, ~] = size(updatedTraceback);
%                 updatedTraceback = sortByRow(updatedTraceback, d1, dist2Ant);
%                 changed = true;
%                 break
%             end
%         end
% 
%     end
% end


updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);
% 
for i=1:20
    [~, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, antPole, postPole, 2, scale_xyz, dist2Post);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
    updatedTraceback = checkRows(updatedTraceback, x, y, z, antPole, postPole);
    NPR = sum(updatedTraceback ~= 0, 2);
    deleteRows = find(NPR <= 2);
    updatedTraceback(deleteRows, :)=[];
end

NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = find(NPR <= 4);
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

changed = true;
counter2 = 0;



while changed && counter2 < 1000
    changed = false;
    counter2 = counter2 + 1;
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);

    [d1, ~] = size(updatedTraceback);
    NPR = sum(updatedTraceback ~= 0, 2);
    heads = updatedTraceback(:, 1);
    tails = zeros(d1, 1);
    for i = 1:d1
        tails(i) = updatedTraceback(i, NPR(i));
    end

    L = 2 * length(tails);
    candidates = zeros(L, 1);
    candidates(1:length(tails)) = tails;
    candidates(length(tails)+1:2*length(tails)) = heads;
    distances = dist2All(candidates, candidates);
    [startPt, endPt] = findClosest(candidates, distances, heads, tails);
    for i = 1:length(startPt)
        if (ismember(startPt(i), heads) && ismember(endPt(i), tails))
            p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
            if true % angleWithPoles < 50/180*pi 
                k1 = find(heads == startPt(i));
                k2 = find(tails == endPt(i));
                row1 = updatedTraceback(k1, 1:NPR(k1));
                row2 = updatedTraceback(k2, 1:NPR(k2));
                zRange1 = min(max(z(row1)) - min(z(row1)), max(z(row2)) - min(z(row2)));
                zRange2 = max(min(max(z(row1)), max(z(row2))) - max(min(z(row1)), min(z(row2))), 0);
                row = zeros(NPR(k1)+NPR(k2), 1);
                row(1:NPR(k1)) = row1;
                row(NPR(k1)+1:end) = row2;
                row = sortByRow(row', 1, dist2Ant);
                if k1 ~= k2 && checkOneRow(row, x, y, z)%&& (zRange2/(zRange1+0.00001) < 0.5)
                    updatedTraceback = connect1Line(updatedTraceback, k1, k2, NPR(k1), NPR(k2), dist2Ant, dist2All);
                    changed = true;
                    break
                end
            end
        elseif (ismember(startPt(i), tails) && ismember(endPt(i), heads))
            p22p1 = [x(startPt(i)), y(startPt(i)), z(startPt(i))] - [x(endPt(i)), y(endPt(i)), z(endPt(i))];
            angleWithPoles = abs(acos(dot(ant2post, p22p1)/(norm(ant2post)*norm(p22p1))));
            if true % angleWithPoles < 50/180*pi 
                k1 = find(tails == startPt(i));
                k2 = find(heads == endPt(i));
                row1 = updatedTraceback(k1, 1:NPR(k1));
                row2 = updatedTraceback(k2, 1:NPR(k2));
                mean1 = [mean(x(row1)), mean(y(row1)), mean(z(row1))];
                mean2 = [mean(x(row2)), mean(y(row2)), mean(z(row2))];
                head1 = updatedTraceback(k1, 1);
                head1 = [mean(x(head1)), mean(y(head1)), mean(z(head1))];
                head2 = updatedTraceback(k2, 1);
                head2 = [mean(x(head2)), mean(y(head2)), mean(z(head2))];
                tail1 = updatedTraceback(k1, NPR(k1));
                tail1 = [mean(x(tail1)), mean(y(tail1)), mean(z(tail1))];
                tail2 = updatedTraceback(k2, NPR(k2));
                tail2 = [mean(x(tail2)), mean(y(tail2)), mean(z(tail2))];  
                if vecnorm(head1-tail2, 2) < vecnorm(head2-tail1, 2)
                    vec1 = head1-tail2;
                else
                    vec1 = head2-tail1;   
                end
                d1 = vecnorm(head1-head2, 2);
                d2 = vecnorm(tail1-tail2, 2);
                d = min(d1, d2);
                row = zeros(NPR(k1)+NPR(k2), 1);
                row(1:NPR(k1)) = row1;
                row(NPR(k1)+1:end) = row2;
                row = sortByRow(row', 1, dist2Ant);
                if k1 ~= k2 && checkOneRow(row, x, y, z) % && vecnorm(vec1, 2) < d %&& (zRange2/(zRange1+0.00001) < 0.5)
                    updatedTraceback = connect1Line(updatedTraceback, k2, k1, NPR(k2), NPR(k1), dist2Ant, dist2All);
                    changed = true;
                    break
                end
            end
        end
    end
end



updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

updatedTraceback = checkRows2(updatedTraceback, x, y, z, antPole, postPole);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

shp = alphaShape(0.8*x, 0.8*y, z);
shp.Alpha = 1.25*shp.Alpha;
[d1, ~] = size(updatedTraceback);
NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = [];
for k1=1:d1
    row1 = updatedTraceback(k1, 1:NPR(k1));
    tf = inShape(shp, x(row1), y(row1), z(row1));
    if sum(tf) >= length(tf)-1
        deleteRows(end+1) = k1;
    end
end
updatedTraceback(deleteRows, :)=[];


[d1, ~] = size(updatedTraceback);
NPR = sum(updatedTraceback ~= 0, 2);
rowLengths = zeros(d1, 1);
for i=1:d1
    L = 0;
    for j =1:NPR(i)-1
        bb1 = [x(updatedTraceback(i, j)), y(updatedTraceback(i, j)), z(updatedTraceback(i, j))];
        bb2 = [x(updatedTraceback(i, j+1)), y(updatedTraceback(i, j+1)), z(updatedTraceback(i, j+1))];
        L = L + vecnorm(bb1-bb2, 2);
    end
    rowLengths(i) = L;
end
% deleteRows = find(rowLengths <= minRowLength);
% updatedTraceback(deleteRows, :)=[];
% withLabel = intersect(1:num_BBs, updatedTraceback);
% withoutLabel = setdiff(1:num_BBs, withLabel);
% 
% updatedTraceback = checkRows(updatedTraceback, x, y, z, antPole, postPole);
% NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = find(NPR <= minBBsInRow2);
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

end