% New version
function [updatedTraceback, withLabel, withoutLabel] = getCiliaryRows3(x, y, z, antPole, postPole, minBBsInRow, minRowLength, scale_xyz)
num_BBs = length(x);
all_pts = [x,y,z].*scale_xyz;
dist2Ant = vecnorm(all_pts - antPole.*scale_xyz,2,2);
dist2Post = vecnorm(all_pts - postPole.*scale_xyz,2,2);

% dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.8, 0.2, 0, 0.0]);
dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.7, 0.3, 0, 0.0]);
[startPt, endPt] = findLink(x, y, z, dist2All, antPole, postPole, dist2Ant, dist2Post);
% matrix showing ciliary rows
OriginalBBsMatrix = link(startPt, endPt);
[d1, ~] = size(OriginalBBsMatrix);
for k=1:d1
    OriginalBBsMatrix = sortByRow(OriginalBBsMatrix, k, dist2Ant);
end
updatedTraceback = OriginalBBsMatrix;

% NPR = sum(updatedTraceback ~= 0, 2);
% deleteRows = find(NPR <= 2);
% updatedTraceback(deleteRows, :)=[];
% withLabel = intersect(1:num_BBs, updatedTraceback);
% withoutLabel = setdiff(1:num_BBs, withLabel);


% changed = true;
ant2post = postPole - antPole;


% 
% updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);


for i=1:1
    [~, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, antPole, postPole, 2, scale_xyz, dist2Post);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);

    [~, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, antPole, postPole, 2, scale_xyz, dist2Post);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
    
    [~, updatedTraceback] = connectLines(updatedTraceback, withLabel, dist2Ant, dist2Post, dist2All, x, y, z, antPole, postPole);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
    
    updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
    
    updatedTraceback = checkRows(updatedTraceback, x, y, z, antPole, postPole);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);

[d1, ~] = size(updatedTraceback);
NPR = sum(updatedTraceback ~= 0, 2);
rowLengths = zeros(d1, 1);
meanH = zeros(d1, 1);
for i=1:d1
    L = 0;
    for j =1:NPR(i)-1
        bb1 = [x(updatedTraceback(i, j)), y(updatedTraceback(i, j)), z(updatedTraceback(i, j))];
        bb2 = [x(updatedTraceback(i, j+1)), y(updatedTraceback(i, j+1)), z(updatedTraceback(i, j+1))];
        L = L + vecnorm(bb1-bb2, 2);
    end
    rowLengths(i) = L;
    meanH(i) = mean(z(updatedTraceback(i, 1:NPR(i))));
end
deleteRows = find(rowLengths <= 10 & meanH > 0.85 * antPole(3));
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);


[~, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, antPole, postPole, 2, scale_xyz, dist2Post);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

[~, updatedTraceback] = connectLines(updatedTraceback, withLabel, dist2Ant, dist2Post, dist2All, x, y, z, antPole, postPole);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

updatedTraceback = checkRows2(updatedTraceback, x, y, z, antPole, postPole);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

end
% for i=1:2
%     [~, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, withoutLabel, dist2Ant, dist2All, x, y, z, antPole, postPole, 2, scale_xyz, dist2Post);
%     withLabel = intersect(1:num_BBs, updatedTraceback);
%     withoutLabel = setdiff(1:num_BBs, withLabel);
% 
%     [~, updatedTraceback] = connectLines(updatedTraceback, withLabel, dist2Ant, dist2Post, dist2All, x, y, z, antPole, postPole);
%     withLabel = intersect(1:num_BBs, updatedTraceback);
%     withoutLabel = setdiff(1:num_BBs, withLabel);
% 
% end
% 
% updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All);
% withLabel = intersect(1:num_BBs, updatedTraceback);
% withoutLabel = setdiff(1:num_BBs, withLabel);
% 
% updatedTraceback = checkRows2(updatedTraceback, x, y, z, antPole, postPole);
% withLabel = intersect(1:num_BBs, updatedTraceback);
% withoutLabel = setdiff(1:num_BBs, withLabel);

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
deleteRows = find(rowLengths <= minRowLength);
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);
% 
% updatedTraceback = checkRows(updatedTraceback, x, y, z, antPole, postPole);
NPR = sum(updatedTraceback ~= 0, 2);
deleteRows = find(NPR <= max(2, minBBsInRow));
updatedTraceback(deleteRows, :)=[];
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);

end