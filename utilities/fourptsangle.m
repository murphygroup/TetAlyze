% -------------------------------------------------------------------------
% [Chongming] 08/20/20
% Calculate the dihedral angle of (pt1, antPole, postPole) and (pt2, antPole, postPole)
% -------------------------------------------------------------------------

function angle = fourptsangle(pt1, pt2, antPole, postPole)

ant2post = postPole - antPole;
% pt1 = pt1 .* scale_xyz;
% pt2 = pt2 .* scale_xyz;
n1 = cross(ant2post, pt1-antPole);
n2 = cross(ant2post, pt2-antPole);
angle = acos(dot(n1,n2)/(norm(n1)*norm(n2)+1e-8));
end
