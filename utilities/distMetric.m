function dist_forward = distMetric(x, y, z, scale, antPole, postPole, dist2Ant, dist2Post, weight)
n = length(x);
dist_forward = zeros(n, n);
dist_Eu = zeros(n, n);
dist2plane = zeros(n, n);
anglesP = zeros(n, n);
angles2Plane = zeros(n, n);
p1 = weight(1);
p2 = weight(2);
p3 = weight(3);
p4 = weight(4);
for i = 1:n
    for j = 1:n
        pt = [x(i), y(i), z(i)];
        pt1 = [x(j), y(j), z(j)];
        if  (dist2Ant(i) < 9) && (dist2Ant(j) < 9)
            dist_Eu(i, j) = distance_pts(pt, pt1,[1, 1, 0.75]);
%             dist_Eu(i, j) = distance_pts(pt, pt1,[1, 1, 1]);
        else
            dist_Eu(i, j) = distance_pts(pt, pt1,[1, 1, 0.5]);
%             dist_Eu(i, j) = distance_pts(pt, pt1,[1, 1, 1]);
        end
            
        if i ~= j
            dist2plane(i, j) = distance_pt2plane(pt1, scale, pt, antPole, postPole);
        else
            dist2plane(i, j) = 0;
        end
        
        ant2post = antPole - postPole;
        p22p1 = pt1 - pt;
        angleWithPoles = acos(abs(dot(ant2post, p22p1))/(norm(ant2post)*norm(p22p1)+1e-9));
        angleWithPoles = angleWithPoles*180/pi;
        anglesP(i,j) = angleWithPoles;
        
        angle = fourptsangle(pt, pt1, antPole, postPole);
        angle = angle*180/pi;
        angles2Plane(i, j) = angle;
%         dist_forward(i, j) = dist_Eu(i, j) + ...
%             20 * dist2plane(i, j);
    end
end

a = mean(dist_Eu(:));
b = mean(dist2plane(:));
c = mean(anglesP(:));
d = mean(angles2Plane(:)); % hurt side
 for i = 1:n
    for j = 1:n

            dist_forward(i, j) = p1*dist_Eu(i, j)/a + p2*dist2plane(i, j)/b + ...
            p3*anglesP(i,j)/c + p4*angles2Plane(i,j)/d; 

    end
end
end