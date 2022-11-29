% -------------------------------------------------------------------------
%% [Mehul Bapat] 23 Feb 2020
%   Vectorized the inner for loop to improve performance and merged the
%   dist_pt_to_plane function with this function.

% [Ben] 12/07/17
% Calculates the minimum distance between each potential BB point, and the
% convex hull created from the set of potential BB points. Returns the
% minimum distances calculated. 
% Note that by 'potential BB point', I am referring to the output of
% 'getPotentialBB_plane'. 
% -------------------------------------------------------------------------

%%  
function minDists = dist2ConvexHull(x, y, z, k, scale_xyz)
% 
% x = x*0.125;
% y = y*0.125;
% z = z*0.3;

x = x*scale_xyz(1);
y = y*scale_xyz(2);
z = z*scale_xyz(3);

minDists = zeros(length(x), 1);

for BB = 1:length(x)
    pt = [x(BB), y(BB), z(BB)];
    
    pt1 = [x(k(:,1)), y(k(:, 1)), z(k(:, 1))];
    pt2 = [x(k(:,2)), y(k(:, 2)), z(k(:, 2))];
    pt3 = [x(k(:,3)), y(k(:, 3)), z(k(:, 3))];
    
    v1_temp = pt1 - pt2;
    v2_temp = pt1 - pt3;
    v_temp = cross(v1_temp,v2_temp,2);
    currDists = abs(dot(pt1 - pt, v_temp,2))./vecnorm(v_temp,2,2);
    
    minDists(BB) = min(currDists);

end
end

% Old Non Vectorized Code
%     for i = 1:size(k, 1)
%         pt1 = [x(k(i, 1)), y(k(i, 1)), z(k(i, 1))];
%         pt2 = [x(k(i, 2)), y(k(i, 2)), z(k(i, 2))];
%         pt3 = [x(k(i, 3)), y(k(i, 3)), z(k(i, 3))];
%         d = dist_pt_to_plane(pt1, pt2, pt3, pt);
%         currDists(i) = d;
%     end
%     