% -------------------------------------------------------------------------
% [Paul] 09/21/20
% Stats of aligned BBs.
% -------------------------------------------------------------------------

function vec = stats (updatedTraceback, x, y, z, antPole, postPole, resultPath, imageID)
num_BBs = length(x);
withLabel = intersect(1:num_BBs, updatedTraceback);
withoutLabel = setdiff(1:num_BBs, withLabel);
numBBs_row = sum(updatedTraceback ~= 0, 2);
[d1, d2] = size(updatedTraceback);


figure(3)
b=bar(numBBs_row, 'yellow');
title('Number of BBs Per Row')
hold off;
saveas(b, join([resultPath, imageID, 'Number of BBs per Row.png']));
writematrix(numBBs_row, join([resultPath, imageID, 'Number of BBs per Row.csv']));

zRange = max(z) - min(z);
zMax = max(z);
zMin = min(z);
intervalDists = zeros(d1, 3);
rowLengths = zeros(d1, 1);
for i = 1:d1
    temp1 = [];
    temp2 = [];
    temp3 = [];
    L = 0;
    for j =1:numBBs_row(i)-1
        bb1 = [x(updatedTraceback(i, j)), y(updatedTraceback(i, j)), z(updatedTraceback(i, j))];
        bb2 = [x(updatedTraceback(i, j+1)), y(updatedTraceback(i, j+1)), z(updatedTraceback(i, j+1))];
        if z(updatedTraceback(i, j)) > zMax - 0.25*zRange && z(updatedTraceback(i, j+1)) > zMax - 0.25*zRange
            temp1(end+1) = vecnorm(bb1 - bb2, 2);
            L = L + temp1(end);
        elseif z(updatedTraceback(i, j)) < zMin + 0.25*zRange && z(updatedTraceback(i, j+1)) < zMin + 0.25*zRange
            temp2(end+1) = vecnorm(bb1 - bb2, 2);
            L = L + temp2(end);
        else
            temp3(end+1) = vecnorm(bb1 - bb2, 2);
            L = L + temp3(end);
        end
    end
    rowLengths(i) = L;
    intervalDists(i, 1) = mean(temp1(:));
    intervalDists(i, 2) = mean(temp2(:));
    intervalDists(i, 3) = mean(temp3(:));
end

figure(4)
subplot(1,2,1)
b=bar(rowLengths, 'cyan');
xlabel('Row Index')
ylabel('Length')
title('Lengths of BB rows (um)')
writematrix(rowLengths, join([resultPath, imageID, 'Lengths of BB rows.csv']));
hold off;

subplot(1,2,2)
h = histogram(rowLengths, 20, 'Normalization', 'probability');
xlabel('Row Length')
ylabel('Frequency')
title('Distribution of BB row lengths')
saveas(b, join([resultPath, imageID, 'Lengths of BB rows.png']));



f = figure(6);
sz = 40;
for i = 1:3
    scatter(1:d1, intervalDists(:, 1) ,sz,'MarkerEdgeColor',[0 .7 .7], 'MarkerFaceColor',[0 .7 .7],  'LineWidth',1.5)
    hold on;
    scatter(1:d1, intervalDists(:, 2) ,sz,'MarkerEdgeColor',[.7 .7 0], 'MarkerFaceColor',[.7 .7 0],  'LineWidth',1.5)
    hold on;
    scatter(1:d1, intervalDists(:, 3) ,sz,'MarkerEdgeColor',[.7 0 .7], 'MarkerFaceColor',[.7 0 .7],  'LineWidth',1.5)
    hold on;
end
avgintervalDists = mean(intervalDists, 'omitnan');
writematrix(intervalDists, join([resultPath, imageID, 'Pairwise BB Distance Raw Data.csv']));


yl = yline(avgintervalDists(1),'--','Average 1','LineWidth',2);
yl.Color = [0 .7 .7];
hold on;
y2 = yline(avgintervalDists(2),'--','Average 2','LineWidth',2);
y2.Color = [.7 .7 0];
hold on;
y3 = yline(avgintervalDists(3),'--','Average 3','LineWidth',2);
y3.Color = [.7 0 .7];
hold on;
legend('Near Anterior Pole','Near Posterior Pole', 'Middle Region');
title('Average distance between each pair of consecutive BBs for each row')
xlabel('Row Index')
ylabel('Distance')
hold off;
saveas(f, join([resultPath, imageID, 'Pairwise BB Distance Scatter Plot.png']));


centroids1 = zeros(d1, 2);
centroids2 = zeros(d1, 2);
centroids3 = zeros(d1, 2);

for i = 1:d1
    temp1 = [];
    temp2 = [];
    temp3 = [];
    for j =1:numBBs_row(i)
        bb1 = [x(updatedTraceback(i, j)), y(updatedTraceback(i, j)), z(updatedTraceback(i, j))];
        if z(updatedTraceback(i, j)) > zMax - 0.25*zRange
            temp1 = cat(3, temp1, bb1(1:2));
        elseif z(updatedTraceback(i, j)) < zMin + 0.25*zRange 
            temp2 = cat(3, temp2, bb1(1:2));
        else
            temp3 = cat(3, temp3, bb1(1:2));
        end  
    end
    if ~isempty(temp1)
    centroids1(i, :) = mean(temp1, 3);
    else
        centroids1(i, 1) = NaN;
        centroids1(i, 2) = NaN;
    end
    if ~isempty(temp2)
    centroids2(i, :) = mean(temp2, 3);
    else
        centroids2(i, 1) = NaN;
        centroids2(i, 2) = NaN;
    end
    if ~isempty(temp3)
    centroids3(i, :) = mean(temp3, 3);
    else
        centroids3(i, 1) = NaN;
        centroids3(i, 2) = NaN;
    end

end

D = pdist(centroids1);
rowDists1 = squareform(D);
for i=1:length(rowDists1)
    rowDists1(i, i) = NaN;
end
rowDists1 = min(rowDists1,[],2);

D = pdist(centroids2);
rowDists2 = squareform(D);
for i=1:length(rowDists2)
    rowDists2(i, i) = NaN;
end
rowDists2 = min(rowDists2,[],2);

D = pdist(centroids3);
rowDists3 = squareform(D);
for i=1:length(rowDists3)
    rowDists3(i, i) = NaN;
end
rowDists3 = min(rowDists3,[],2);

f = figure(7);
sz = 40;
for i = 1:3
    scatter(1:d1, rowDists1 ,sz,'MarkerEdgeColor',[0 .7 .7], 'MarkerFaceColor',[0 .7 .7],  'LineWidth',1.5)
    hold on;
    scatter(1:d1, rowDists2 ,sz,'MarkerEdgeColor',[.7 .7 0], 'MarkerFaceColor',[.7 .7 0],  'LineWidth',1.5)
    hold on;
    scatter(1:d1, rowDists3 ,sz,'MarkerEdgeColor',[.7 0 .7], 'MarkerFaceColor',[.7 0 .7],  'LineWidth',1.5)
    hold on;
end


yl = yline(mean(rowDists1, 'omitnan'),'--','Average 1','LineWidth',2);
yl.Color = [0 .7 .7];
hold on;
y2 = yline(mean(rowDists2, 'omitnan'),'--','Average 2','LineWidth',2);
y2.Color = [.7 .7 0];
hold on;
y3 = yline(mean(rowDists3, 'omitnan'),'--','Average 3','LineWidth',2);
y3.Color = [.7 0 .7];
hold on;
legend('Near Anterior Pole','Near Posterior Pole', 'Middle Region');
title('Distance between Clliary BB Rows')
xlabel('Row Index')
ylabel('Distance')
hold off;
saveas(f, join([resultPath, imageID, 'BB Row Spacing Scatter Plot.png']));

tempU = zeros(length(rowDists1), 3);
tempU(:, 1) = rowDists1;
tempU(:, 2) = rowDists2;
tempU(:, 3) = rowDists3;

writematrix(tempU, join([resultPath, imageID, 'BB Row Spacing Raw Data.csv']));

vec = [avgintervalDists(1), avgintervalDists(2), avgintervalDists(3), mean(rowDists1, 'omitnan'), mean(rowDists2, 'omitnan'), mean(rowDists3, 'omitnan'), ...
    mean(rowLengths), std(rowLengths)];
end