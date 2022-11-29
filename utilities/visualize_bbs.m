% -------------------------------------------------------------------------
% [Chongming] 08/20/20
% 3D/2D visualization of basal bodies of a cell for given a given points
% matrix
% -------------------------------------------------------------------------

function visualize_bbs(x, y, z, oa_x, oa_y, oa_z, pts, asignedPoints, unsignedPoints, antPole, postPole, scale_xyz)

% [x, y, z] = newCoor(x, y, z, antPole, postPole);
x=x*scale_xyz(1);
y=y*scale_xyz(2);
z=z*scale_xyz(3);
oa_x=oa_x*scale_xyz(1);
oa_y=oa_y*scale_xyz(2);
oa_z=oa_z*scale_xyz(3);
antPole = antPole.*scale_xyz;
postPole = postPole.*scale_xyz;

% K = convhull(x(asignedPoints), y(asignedPoints), z(asignedPoints));
shp = alphaShape(0.9*x(asignedPoints), 0.9*y(asignedPoints), 0.98*z(asignedPoints));
shp.Alpha = 1.5 * shp.Alpha;
% plot mesh of convex hull of BB positions
% trimesh(K, 0.8*x(asignedPoints), 0.8*y(asignedPoints), 0.97*z(asignedPoints), 'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.8, ...
%     'EdgeColor', [150 150 150]/255);
plot(shp,'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.8,  'EdgeColor', [150 150 150]/255, 'DisplayName', 'Cell Body');
hold on
colors = {'green', 'magenta', 'cyan', 'blue', 'yellow', 'red', 'black'};
lineStyles = {'-', '--'};
markers = {'o', '*', 'x', 's', 'd', 'p', 'h'};
[num_rows, ~] = size(pts);
xy_centroids = zeros(num_rows, 1);
for i = 1:num_rows
    row = pts(i, :);
    row = row(row ~= 0);
    xy_centroids(i) = atan2(mean(y(row)), mean(x(row)));
end
[~, color_order] = sort(xy_centroids, 'ascend');
for i = 1:num_rows
    row = pts(i, :);
    row = row(row ~= 0);
    lineStyle = lineStyles{fix(i/49) + 1};
    temp = mod(i, 49);
    color = colors{mod(temp, 7) + 1};
    marker = markers{fix(temp/7) + 1};
    plot3(x(row), y(row), z(row),'MarkerSize',7, 'Color', color, 'LineStyle', lineStyle, 'Marker', marker, 'LineWidth', 1.5, 'DisplayName',strcat('Row ',int2str(i)));
    hold on
    axis equal;
    
end
scatter3(x(unsignedPoints), y(unsignedPoints), z(unsignedPoints),13, '+', 'MarkerEdgeColor','red', 'DisplayName', 'Unassigned')


scatter3(oa_x, oa_y, oa_z, 10, 'MarkerEdgeColor','yellow', 'MarkerFaceColor', 'yellow',  'DisplayName', 'OA')
hold on
axis equal;

scatter3(antPole(1), antPole(2), antPole(3), 15, 'd','DisplayName', 'Anterior Pole');
scatter3(postPole(1), postPole(2), postPole(3), 15,  'DisplayName', 'Posterior Pole');
legend
% dim = [.2 .6 .3 .3];
% str ='Red point refers to the unassigned basel body, and yellow point refers to the oral apparatus';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
end
