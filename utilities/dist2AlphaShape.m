function minDists = dist2AlphaShape(x, y, z, k, scale_xyz)

d1 = length(x);
d2 = length(k);
distMtx = zeros(d1, d2);
for i=1:d1
    pt = [x(i), y(i), z(i)];
    for j=1:d2
        pt1 = k(j,:);
        x1=pt(1);
        y1=pt(2);
        z1=pt(3);

        x2=pt1(1);
        y2=pt1(2);
        z2=pt1(3);
        distMtx(i, j) = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
    end
end

minDists = min(distMtx,[],2);
end