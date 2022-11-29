function updatedTraceback = checkRows2(updatedTraceback, x, y, z, antPole, postPole)

changed = true;
while changed
    changed = false; 
    [d1, ~] = size(updatedTraceback);
    NPR = sum(updatedTraceback ~= 0, 2);
    for k=1:d1
        if NPR(k) > 2
            for i=2:NPR(k)-1
                p1 = [x(updatedTraceback(k, i-1)), y(updatedTraceback(k, i-1)), z(updatedTraceback(k, i-1))];
                p2 = [x(updatedTraceback(k, i)), y(updatedTraceback(k, i)), z(updatedTraceback(k, i))];
                p3 = [x(updatedTraceback(k, i+1)), y(updatedTraceback(k, i+1)), z(updatedTraceback(k, i+1))];
                vec1 = p2 - p1;
                vec2 = p2 - p3;
                angle = acos(dot(vec1, vec2)/(norm(vec1)*norm(vec2)));
                angle = angle * 180 /pi;
                if angle < 40
                    updatedTraceback(k,i:NPR(k)-1) = updatedTraceback(k, i+1:NPR(k));
                    updatedTraceback(k, NPR(k)) = 0;
                    changed = true;
                    break
                end
            end
        end
        if changed
            break    
        end 
    end
end

ant2post = antPole - postPole;
changed = true;
while changed
    changed = false; 
    [d1, ~] = size(updatedTraceback);
    NPR = sum(updatedTraceback ~= 0, 2);
    for k=1:d1
        meanDist = 0;
        stdDist = 0;
        temp = [];
        for j =1:NPR(k)-1
            bb1 = [x(updatedTraceback(k, j)), y(updatedTraceback(k, j)), z(updatedTraceback(k, j))];
            bb2 = [x(updatedTraceback(k, j+1)), y(updatedTraceback(k, j+1)), z(updatedTraceback(k, j+1))];
            temp(end+1) = vecnorm(bb1 - bb2, 2);     
        end
        meanDist = mean(temp(:));
        stdDist = std(temp(:));
        if NPR(k) > 2
            for i=1:NPR(k)-1
                p1 = [x(updatedTraceback(k, i)), y(updatedTraceback(k, i)), z(updatedTraceback(k, i))];
                p2 = [x(updatedTraceback(k, i+1)), y(updatedTraceback(k, i+1)), z(updatedTraceback(k, i+1))];
                vec1 = p2 - p1;
                dist = vecnorm(vec1, 2);
                angleLocal = acos(dot(vec1, ant2post)/(norm(vec1)*norm(ant2post)));
                angleLocal = angleLocal/pi*180;
                if (dist > 3*meanDist || dist > 8) && max(z(updatedTraceback(k, i+1)), z(updatedTraceback(k, i))) > 0.75 * antPole(3)
                    if i == 1
                        updatedTraceback(k,1:NPR(k)-1) = updatedTraceback(k, 2:NPR(k));
                        updatedTraceback(k,NPR(k)) = 0;
                    elseif i == NPR(k) -1
                        updatedTraceback(k, NPR(k)) = 0;
                    else
                        updatedTraceback(end+1,1:NPR(k) - i) = updatedTraceback(k, i+1:NPR(k));
                        updatedTraceback(k,i+1:NPR(k)) = 0;
                    end
                    changed = true;
                    break
                end
            end
        end
        if changed
            break
        end
    end   
end

end