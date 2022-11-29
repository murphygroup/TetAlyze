function D = dist2line(x,y,z)
%dist2line Returns the distance matrix for each point
%% [Mehul Bapat] 26-Feb-2020
%   Calculates the distance between a point and the line passing through
%   the mean of all data and the mean of all data but at a distance of
%   minimum z.

xm = mean(x);
ym = mean(y);
zm = mean(z);

pt1 = [xm,ym,zm]; %mean of all data
pt2 = [xm,ym,min(z)]; %in the first plane of data
pt = [x(:),y(:),z(:)];

v1 = pt2-pt1;
v2 = pt-pt1;

D = abs(vecnorm(v2,2,2).*sin( acos(  dot(v1.*ones(size(v2)),v2,2)./(vecnorm(v1).*vecnorm(v2,2,2))  ) ));

end