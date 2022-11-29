% -------------------------------------------------------------------------
% [Ben] 3/28/18
% Takes in the x,y,z coordinates of cortical BBs and returns estimated
% anterior and posterior poles of the cell.
% -------------------------------------------------------------------------

function [pole1, pole2] = findPoles(x,y,z,scale_xyz)
x=x*scale_xyz(1);
y=y*scale_xyz(2);
z=z*scale_xyz(3);
pts = [x y z];

numBBs=length(x);

BBsDistanceMatrix = zeros(numBBs);
maxDistanceForEachBB = zeros(numBBs,1);
for i=1:numBBs
    tempdist = vecnorm(pts-pts(i,:),2,2);
    BBsDistanceMatrix(i,:) = tempdist.';
    maxDistanceForEachBB(i,1)=max(BBsDistanceMatrix(i,:));
end

% use numBBFindPoles of max distance to determine poles
numBBFindPoles = 20;
maxDistanceForEachBB_sorted=sort(maxDistanceForEachBB,'descend');

BBPairWithMaxDistance_i=zeros(numBBFindPoles,1);
BBPairWithMaxDistance_j=zeros(numBBFindPoles,1);

% pg. 18 of Jingyi's paper. Get pairwise distance matrix for all detected
% cortical BBs. Find the 10 cortical BB pairs that produce the greatest
% inter-BB distance. 
for k=1:numBBFindPoles
    for i=1:numBBs
        for j=1:numBBs
            if BBsDistanceMatrix(i,j)==maxDistanceForEachBB_sorted(k)
                BBPairWithMaxDistance_i(k)=i;
                BBPairWithMaxDistance_j(k)=j;
            end
        end
    end
end

BBPairWithMaxDistance=[BBPairWithMaxDistance_i;BBPairWithMaxDistance_j];
% extracts unique values in matrix
BBPairWithMaxDistance=unique(BBPairWithMaxDistance);

idxBBPairWithMaxDistance=zeros(length(BBPairWithMaxDistance),3);
for i=1:length(BBPairWithMaxDistance)
    idxBBPairWithMaxDistance(i,1)=x(BBPairWithMaxDistance(i));
    idxBBPairWithMaxDistance(i,2)=y(BBPairWithMaxDistance(i));
    idxBBPairWithMaxDistance(i,3)=z(BBPairWithMaxDistance(i));
end

% pg. 18 of Jingyi's paper. Apply k-means to basal bodies forming the 10
% greatest inter-BB distances. Anterior and posterior pole defined as the
% centroids of two clusters. 
numOfGroup=2;
[~, groupCenter]=kmeans(idxBBPairWithMaxDistance,numOfGroup);

% -------------------------------------------------------------------------
% [Chongming] 8/20/20
% The pole show be divided by scale_xyz
% -------------------------------------------------------------------------

% pole1 = groupCenter(1, :) ./ [0.125 0.125 0.3];
% pole2 = groupCenter(2, :) ./ [0.125 0.125 0.3];

pole1 = groupCenter(1, :) ./ scale_xyz;
pole2 = groupCenter(2, :) ./ scale_xyz;

end