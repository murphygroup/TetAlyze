function sig = checkOneRow(row, x, y, z)
sig = true;
tol = max(0.1*length(row),2);
n = 0;
if length(row) > 2
    for i=2:length(row)-1
        p1 = [x(row(i-1)), y(row(i-1)), z(row(i-1))];
        p2 = [x(row(i)), y(row(i)), z(row(i))];
        p3 = [x(row(i+1)), y(row(i+1)), z(row(i+1))];
        vec1 = p2 - p1;
        vec2 = p2 - p3;
        angle = acos(dot(vec1, vec2)/(norm(vec1)*norm(vec2)));
        angle = angle * 180 /pi;
        d1 = vecnorm(vec1, 2);
        d2 = vecnorm(vec2, 2);
        
        if angle < 70 %&& min(d1, d2) > 0.36
            n = n + 1;
            if n >= tol
                sig = false;
                break
            end
        end
    end
end
end