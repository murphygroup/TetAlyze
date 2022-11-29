function [pts] = getpts3d(x, y, z, n)
    %% getpts3d Select points from a 3D scatter plot by clicking on plot
    % x - Vector of X coordinates
    % y - Vector of Y coordinates
    % z - Vector of Z coordinates
    % n - Number of points needed to be selected
    % pts - Returns a n x 3 matrix of selected points

    h = figure;
    scatter3(x,y,z);
    pts = zeros(n,3);
    datacursormode on
    dcm_obj = datacursormode(h);

    for i=1:n
        display(sprintf('Click on figure for point %d',i))
        waitforbuttonpress;
        f = getCursorInfo(dcm_obj);
        pts(i,:) = f.Position;
    end

end