function immBBnum = intensityAnalysis(cort_x, cort_y, cort_z, BW, I, firstRow, resultPath)
    x_ = 2*(cort_x - min(cort_x))/(max(cort_x) - min(cort_x))-1;

    y_ = 2*(cort_y - min(cort_y))/(max(cort_y) - min(cort_y)) - 1;
    z_ = 2*(cort_z - min(cort_z))/(max(cort_z) - min(cort_z)) - 1;
    theta = atan(y_./x_);
    for i=1:length(x_)
        if x_(i) < 0
            theta(i) = theta(i) + pi;
        end
    end

    [theta_,rho,z] = cart2pol(x_,y_,z_);
    r = sqrt(x_.^2+y_.^2+z_.^2);
    phi = acos(z_./r);    
    values = zeros(length(cort_x), 1);
    coordinates = [cort_x, cort_y, cort_z];
    colors = zeros(length(cort_x), 1);
    theta = theta - mean(theta(nonzeros(firstRow)));
    BW = cast(BW, 'uint32');
    index = [];
    for i=1:length(cort_x)
        Idx = knnsearch(coordinates, [cort_x(i), cort_y(i), cort_z(i)], 'K', 11);
%         int = zeros(10, 1);
%         for j=2:11
%             int(j-1) = I(BW(Idx(j), 2), BW(Idx(j), 1), BW(Idx(j), 3));
%         end
%         [h,p] = ztest(I(BW(i, 2), BW(i, 1), BW(i, 3)), mean(int), std(int), 'Tail','left');
        r_ = I(BW(Idx(2), 2), BW(Idx(2), 1), BW(Idx(2), 3)) / I(BW(i, 2), BW(i, 1), BW(i, 3));
        colors(i) = r_ > 2;
        if r_ > 2
            index(end+1) = i;
        end
        values(i) = r_;
    end
            
    writematrix([theta, phi, r, values], join([resultPath, '/', 'intensityAnalysis.csv']));
    figure(8);
    scatter3(cos(theta).*sin(phi), sin(theta).*sin(phi), cos(phi), 20, colors);

    colormap(jet);
    colorbar;

    idx = find(colors);
    antIdx = find(cos(phi)>0.5);
    nAnt = length(nonzeros(ismember(idx, antIdx)));

    postIdx = find(cos(phi)<-0.5); 
    nPost = length(nonzeros(ismember(idx, postIdx)));
    % 
    medIdx1 = find(-0.5<cos(phi) & cos(phi)<0); 
    nMed1 = length(nonzeros(ismember(idx, medIdx1)));

    medIdx2 = find(0<cos(phi) & cos(phi)<0.5); 
    nMed2 = length(nonzeros(ismember(idx, medIdx2)));
    
    newBBSliceSave(BW(index, :), I, resultPath, '/')
    
    immBBnum = [nPost, nMed1, nMed2, nAnt];


        
end