% -------------------------------------------------------------------------
% [Ben] 12/19/17 (adapted from Jingyi's code, by Ben)
% Returns the optimized ciliary row configuration, as well as the
% corresponding F-score
% -------------------------------------------------------------------------

function [bestAlignment, withLabel, withoutLabel] = getCiliaryRows(x, y, z, oa_x, oa_y, oa_z, antPole, postPole, minBBsInRow,scale_xyz)
% To identify the poles of the cell, the maximum distance between each BB
% and every other BB is determined and a list of the 10 greatest distances
% is created. These distances correspond to the 10 most anterior BBs and
% the 10 most posterior BBs and the centroid of each of these clusters of
% BBs is taken to be the anterior and posterior pole respectively.

num_BBs = length(x);
all_pts = [x,y,z].*scale_xyz;
dist2Ant = vecnorm(all_pts - antPole.*scale_xyz,2,2);
dist2Post = vecnorm(all_pts - postPole.*scale_xyz,2,2);
% 
% BBs are assigned a position within a cortical row or kinety by using a
% metric that minimizes the distance between each BB and its partner while
% also minimizing the distance between its partner and a plane comprising
% the 3D coordinates of the BB maxima, the anterior pole and the posterior
% pole.
dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.4,0.5,0,0.1]);
[startPt, endPt] = findLink(x, y, z, dist2All, antPole, postPole, dist2Ant, dist2Post);
% matrix showing ciliary rows
OriginalBBsMatrix = link(startPt, endPt);


% -------------------------------------------------------------------------

% withLabel contains BBs that belong to a particular ciliary row
% withoutLabel contains BBs that were not assigned to any row
withLabel = intersect(1:num_BBs, OriginalBBsMatrix);
withoutLabel = setdiff(1:num_BBs, withLabel);

updatedTraceback = OriginalBBsMatrix;
counter = 1;
% visualize_bbs(x, y, z, 3*[], 3*[],[], ...
%     updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
bestScore = -1;

numBBs = length(withLabel) + length(withoutLabel);
label = zeros(numBBs, 1);
numBBs_row = sum(updatedTraceback ~= 0, 2);
[numRows, ~] = size(updatedTraceback);
for i = 1:numRows
    bbs = updatedTraceback(i, 1:numBBs_row(i));
    label(bbs) = i;
end
n_rounds = 1000;
dist2AllForConnection = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.7,0.2,0,0.1]);
dist2All = distMetric(x, y, z, scale_xyz, antPole, postPole, dist2Ant, dist2Post, [0.6,0.3,0,0.1]);
while counter < n_rounds   
    for k=1:5
        [label, updatedTraceback] = getIniLabel(updatedTraceback, withLabel, ...
                                            withoutLabel, dist2Ant, dist2All, x, y, z, ...
                                            antPole, postPole, minBBsInRow,scale_xyz,dist2Post);
        withLabel = intersect(1:num_BBs, updatedTraceback);
        withoutLabel = setdiff(1:num_BBs, withLabel);
    end
%     visualize_bbs(x, y, z, 3*[], 3*[],[], ...
%         updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
    
    [d, ~] = size(updatedTraceback);
    % swap rows
    for i=1:100
       r = randi([1,d],2,1);
       updatedTraceback(r,:) = updatedTraceback(r([2,1]),:);
    end
    
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
    
    [label, updatedTraceback] = connectLines(updatedTraceback, label, dist2Ant, dist2Post, dist2AllForConnection, x, y, z, antPole, postPole);
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
%     visualize_bbs(x, y, z, 3*[], 3*[],[], ...
%         updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);
    [d, ~] = size(updatedTraceback);
    if length(withLabel)/d > bestScore && counter >= 5
        bestAlignment = updatedTraceback;
        bestScore = length(withLabel)/d;
    end
    
    if counter ~= n_rounds - 1
        numBBs_row = sum(updatedTraceback ~= 0, 2);
        deleteRows = find(numBBs_row <= 7);
%         temp = [];
%         for k=1:length(deleteRows)
%             r = randi([0 2]);
%             if r==1
%                 temp(end+1)=k;
%             end
%         end
%         deleteRows(temp) = [];
        updatedTraceback(deleteRows, :)=[];      
        updatedTraceback = deleteAssignedPoints(65, updatedTraceback,scale_xyz, x, y, z, antPole, postPole);
    end
%     visualize_bbs(x, y, z, 3*[], 3*[],[], ...
%         updatedTraceback, withLabel, withoutLabel, antPole, postPole, scale_xyz);

   
    [d, ~] = size(updatedTraceback);
    % swap rows
    for i=1:100
       r = randi([1,d],2,1);
       updatedTraceback(r,:) = updatedTraceback(r([2,1]),:);
    end
    withLabel = intersect(1:num_BBs, updatedTraceback);
    withoutLabel = setdiff(1:num_BBs, withLabel);
    counter = counter + 1;
end

numBBs_row = sum(bestAlignment ~= 0, 2);
deleteRows = find(numBBs_row <= 3);
bestAlignment(deleteRows, :)=[];      
withLabel = intersect(1:num_BBs, bestAlignment);
withoutLabel = setdiff(1:num_BBs, withLabel);                 
end
                                     

