function vec = analysis(updatedTraceback, x, y, z, antPole, postPole, resultPath, imageID, filename, immBBnum)
% filename = 'Summary.txt'
fid = fopen(join([resultPath, imageID, filename]),'a');


nBBs = length(x)
fprintf(fid, '# identified BBs: %d\n', nBBs);

withLabel = intersect(1:nBBs, updatedTraceback);
withoutLabel = setdiff(1:nBBs, withLabel);
nAssignedBBs = length(withLabel)
fprintf(fid, '# assigned BBs: %d\n', nAssignedBBs);

[d, ~] = size(updatedTraceback);
nRows = d
fprintf(fid, '# BB rows: %d\n', nRows);


NPR = sum(updatedTraceback ~= 0, 2);

nBBsPerRow = mean(NPR)
fprintf(fid, 'Avergae number of BBs per row: %.3f\n', nBBsPerRow);


nBBsPerRowSTDEV = std(NPR)
v = stats(updatedTraceback, x, y, z, antPole, postPole, resultPath, imageID);

xRange = max(x(withLabel)) - min(x(withLabel));

yRange = max(y(withLabel)) - min(y(withLabel));
zRange = max(z(withLabel)) - min(z(withLabel))
fprintf(fid, 'Cell height(um): %.3f\n', zRange);




step = 2;
longaxis = -1;


data =[x y];
[coeff,score,latent] = pca(data);
w1 = (max(score(1:end, 1))-min(score(1:end, 1)));
w2 = (max(score(1:end, 2))-min(score(1:end, 2)));
fprintf('Cell widths(um): %.3f, %.3f\n', w1, w2);
fprintf(fid, 'Cell widths(um): %.3f, %.3f\n', w1, w2);
    
ellipsoidVolume = 4/3 * pi * w1 * w2 * zRange / 8;


qualified = ones(length(withLabel),1);
step = 1;
cort_x = x(withLabel);
cort_y = y(withLabel);
cort_z = z(withLabel);
for i = min(cort_z):step:max(cort_z)
    [idx1, ~] = find(i <= cort_z);
    [idx2, ~] = find(cort_z < min(i+ step, max(cort_z)));
    idx = intersect(idx1, idx2);
    shp1 = alphaShape(0.85*cort_x(idx), 0.85*cort_y(idx), cort_z(idx));
    if shp1.Alpha ~= inf && (i >= min(cort_z) + 1.5 && i <= max(cort_z) - 1.5)
        BBs = [0.9*cort_x(idx), 0.9*cort_y(idx), cort_z(idx)];
        if i < min(cort_z)
            BBs(end+1,1:3) = [0, 0, i];
            BBs(end+1,1:3) = [1, 1, i];
            BBs(end+1,1:3) = [1, -1, i];
            BBs(end+1,1:3) = [-1, 1, i];
            BBs(end+1,1:3) = [-1, -1, i];
        end
        if i+ 4 < max(cort_z)
            BBs(end+1,1:3) = [0, 0, i+step];
            BBs(end+1,1:3) = [1, 1, i+step];
            BBs(end+1,1:3) = [1, -1, i+step];
            BBs(end+1,1:3) = [-1, 1, i+step];
            BBs(end+1,1:3) = [-1, -1, i+step];
        end
        shp1 = alphaShape(BBs(:, 1), BBs(:, 2), BBs(:, 3));
        shp1.Alpha = 2 * shp1.Alpha;
        tf = inShape(shp1, cort_x(idx), cort_y(idx), cort_z(idx));
    elseif i < min(cort_z) + 2 || i > max(cort_z) - 2
        tf = zeros(length(idx), 1);
    else
        tf = ones(length(idx), 1);
    end

    for j = 1:length(idx)
        if tf(j) ~= 0 
            qualified(idx(j)) = 0;
        end
    end
end


cort_x = cort_x(find(qualified));
cort_y = cort_y(find(qualified));
cort_z = cort_z(find(qualified));


shp = alphaShape(cort_x, cort_y, cort_z);

shp.Alpha = 1.25*shp.Alpha;
% figure(6);
% plot(shp);

CellVol = volume(shp)
fprintf(fid, 'Cell volume(um^3): %.3f\n', CellVol);

VolumeDeficit = (ellipsoidVolume - CellVol)/ellipsoidVolume
fprintf(fid, 'VolumeDeficit: %.3f\n', VolumeDeficit);
% [X,Y,Z] = ellipsoid(0,0,0,xRange,yRange,zRange);
% figure(7);
% surf(X,Y,Z);
% axis equal
fprintf(fid, 'Average neighbor BB pairwise distance(um) (anterior, medial, posterior): %.3f, %.3f, %.3f\n', v(1), v(3), v(2));
fprintf(fid, 'Average BB row pairwise distance(um) (anterior, medial, posterior): %.3f, %.3f, %.3f\n', v(4), v(6), v(5));
fprintf('Average neighbor BB pairwise distance(um) (anterior, medial, posterior): %.3f, %.3f, %.3f\n', v(1), v(3), v(2));
fprintf('Average BB row pairwise distance(um) (anterior, medial, posterior): %.3f, %.3f, %.3f\n', v(4), v(6), v(5));
fprintf('Mean and standard deviation of ciliary row lengths: %.3f, %.3f\n', v(7), v(8));
fprintf(fid, 'Mean and standard deviation of ciliary row lengths: %.3f, %.3f\n', v(7), v(8));
% fprintf(fid, 'Number of immatured BB in posterior, posterior medial, anterior medial, and anterior regions: %d, %d, %d, %d\n', immBBnum(1), immBBnum(2), immBBnum(3), immBBnum(4));
% fprintf('Number of immatured BB in posterior, posterior medial, anterior medial, and anterior regions: %d, %d, %d, %d\n', immBBnum(1), immBBnum(2), immBBnum(3), immBBnum(4));
fclose(fid);
vec = [nBBs, nAssignedBBs, nRows, nBBsPerRow, nBBsPerRowSTDEV, xRange, yRange, zRange, w1, w2, CellVol, VolumeDeficit,...
    v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8)];
end