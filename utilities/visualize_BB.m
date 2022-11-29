% -------------------------------------------------------------------------
% [Ben] 05/01/18 (Written by Ben)
% 3D/2D visualization of basal bodies of a cell for given .nd2 file in
% /nd2_files. The cortical BBs are colored 'red' and the oral apparatus BBs
% are colored 'green'. The line anterior-posterior axis is drawn as a
% dotted white line. This function takes in a path to a .nd2 file, and 
% takes about 13s to run. Depending on what the parameter 'slice' is set
% to, either a 3D mesh plot (of the convex hull of BB positions) or some 2D
% slice of that 3D plot is produced.
% Have not implemented option to control aspect ratio for data.
% -------------------------------------------------------------------------


function visualize_BB(cort_x, cort_y, cort_z, oa_x, oa_y, oa_z, antPole, postPole, slice)
shp = alphaShape(0.9*cort_x, 0.9*cort_y, 0.98*cort_z);
shp.Alpha = 1.5 * shp.Alpha;
% plot mesh of convex hull of BB positions
% trimesh(K, 0.8*x(asignedPoints), 0.8*y(asignedPoints), 0.97*z(asignedPoints), 'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.8, ...
%     'EdgeColor', [150 150 150]/255);
plot(shp,'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.8,  'EdgeColor', [150 150 150]/255);
hold on
axis equal;
if strcmp(slice, 'XYZ')
%     set(gca,'BoxStyle','full','Box','on');

    % plot cortical BBs as red filled circles 
    scatter3(cort_x, cort_y, cort_z, 160, 'r.');
    hold on
    % plot oral apparatus BBs as green filled circles
    scatter3(oa_x, oa_y, oa_z, 160, 'g.');
    hold on
    % plot yellow straight line between the anterior and posterior poles
    plot3([postPole(1), antPole(1)], [postPole(2), antPole(2)], ...
        [postPole(3), antPole(3)], '-', 'Color', 'yellow', 'LineWidth', 3);
%     % plot OA centroid
%     scatter3(OA(1), OA(2), OA(3), 2000, 'y.');
%     % plot perpendicular from OA centroid to anterior-posterior axis
%     plot3([OA(1) antPole(1)], [OA(2) antPole(2)], [OA(3) OA(3)], ...
%         '-', 'Color', [1, 1, 1], 'LineWidth', 2);

    hold off
%     zlabel(slice(3));
elseif strcmp(slice, 'XY')
    K = convhull(x, y);
    fill(x(K), y(K), [0.5 0.5 1]);
    hold on
    scatter(cort_x, cort_y, 80, 'r.');
    scatter(oa_x, oa_y, 80, 'g.');
    hold off
elseif strcmp(slice, 'XZ')
    K = convhull(x, z);
    fill(x(K), z(K), [0.5 0.5 1]);
    hold on
    scatter(cort_x, cort_z, 80, 'r.');
    scatter(oa_x, oa_z, 80, 'g.');
    daspect([1 (0.125/0.3) 1]);
    hold off
elseif strcmp(slice, 'YZ')
    K = convhull(y, z);
    fill(y(K), z(K), [0.5 0.5 1]);
    hold on
    scatter(cort_y, cort_z, 80, 'r.');
    scatter(oa_y, oa_z, 80, 'g.');
    daspect([1 (0.125/0.3) 1]);
    hold off
end
title(slice);
xlabel(slice(1));
% [Mehul [20 Feb 2020]] Temporarily changed value in ylabel slice from 2 to 1
ylabel(slice(2));
end

