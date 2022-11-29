% Paul 11/12
function updatedTraceback = refinement(updatedTraceback, x, y, z, dist2Ant, antPole, postPole, dist2All)

changed = true;
counter = 0;

ignored = [];
while changed
    if counter > 500
        break
    end
    changed = false;
    counter = counter + 1;
    [d1, ~] = size(updatedTraceback);
    NPR = sum(updatedTraceback ~= 0, 2);
    heads = updatedTraceback(:, 1);
    tails = zeros(d1, 1);
    for i = 1:d1
        tails(i) = updatedTraceback(i, NPR(i));
    end
    distances = zeros(d1, d1);
    closest = zeros(d1, 1);
    for i=1:d1
        for j=1:d1
            if i ~= j
                distances(i, j) = dist2All(heads(i), tails(j));
            else
                distances(i, j) = 1e10;
            end
        end
        [v, idx] = sort(distances(i, :), 'ascend');
        closest(i) = v(1);
    end
    [rows, ] = find(closest < 0.5);
    minDists = closest(rows);
    minDists = sort(minDists, 'ascend');
    for q =1 : length(minDists)
        if minDists(q) < 0.5
            [k1, k2] = find(distances == minDists(q));
            if length(k1) >1
                k1 = k1(1);
            end
            if length(k2) >1
                k2 = k2(1);
            end
            row1 = updatedTraceback(k1, 1:NPR(k1));
            row2 = updatedTraceback(k2, 1:NPR(k2));
%             row = zeros(NPR(k1)+NPR(k2), 1);
%             row(1:NPR(k1)) = row1;
%             row(NPR(k1)+1:end) = row2;
%             row = sortByRow(row', 1, dist2Ant);
            if true % checkOneRow(row, x, y, z)
                mean1 = [mean(x(row1)), mean(y(row1)), mean(z(row1))];
                mean2 = [mean(x(row2)), mean(y(row2)), mean(z(row2))];
                head1 = updatedTraceback(k1, 1);
                head1 = [mean(x(head1)), mean(y(head1)), mean(z(head1))];
                head2 = updatedTraceback(k2, 1);
                head2 = [mean(x(head2)), mean(y(head2)), mean(z(head2))];
                tail1 = updatedTraceback(k1, NPR(k1));
                tail1 = [mean(x(tail1)), mean(y(tail1)), mean(z(tail1))];
                tail2 = updatedTraceback(k2, NPR(k2));
                tail2 = [mean(x(tail2)), mean(y(tail2)), mean(z(tail2))];  
                if vecnorm(head1-tail2, 2) < vecnorm(head2-tail1, 2)
                    vec1 = head1-tail2;
                else
                    vec1 = head2-tail1;   
                end
                vec3 = mean1-mean2;
                d1 = vecnorm(head1-head2, 2);
                d2 = vecnorm(tail1-tail2, 2);
                d = min(d1, d2);
                vec2 = antPole-postPole;
                
                if vecnorm(vec1, 2) < d
                    localDist = dist2All(row1, row2);
                    minDist = min(localDist(:));
                    [p1, p2] = find(localDist == minDist);
                    if p1 > p2
                        continue
                    end
                    BB1 = [x(updatedTraceback(k1, p1)), y(updatedTraceback(k1, p1))];
                    BB2 = [x(updatedTraceback(k2, p2)), y(updatedTraceback(k2, p2))];
                    vec4 = BB1 - BB2;
                    d2 = vecnorm(vec4, 2);
                    angle = acos(abs(dot(vec1, vec2))/(norm(vec1)*norm(vec2)));
                    angle = angle/pi*180;
                    angle2 = fourptsangle(mean1, mean2, antPole, postPole);
                    angle2 = angle2/pi*180;
                    row = zeros(NPR(k1)+p2-p1+1, 1);
                    row(1:NPR(k1)-p1+1) = updatedTraceback(k1, p1:NPR(k1));
                    row(NPR(k1)-p1+2:end) = updatedTraceback(k2, 1:p2);
                    row = sortByRow(row', 1, dist2Ant);
                    
                    sig = false;
                    if p1 < NPR(k1)-2 && p2 > 3
                        BB1 = [x(updatedTraceback(k1, p1)), y(updatedTraceback(k1, p1)), z(updatedTraceback(k1, p1))];
                        BB2 = [x(updatedTraceback(k2, p2)), y(updatedTraceback(k2, p2)), z(updatedTraceback(k2, p2))];
                        BB3 = [mean(x(updatedTraceback(k1, p1+1:p1+3))), mean(y(updatedTraceback(k1, p1+1:p1+3))), mean(z(updatedTraceback(k1, p1+1:p1+3)))];
                        BB4 = [mean(x(updatedTraceback(k2, p2-3:p2-1))), mean(y(updatedTraceback(k2, p2-3:p2-1))), mean(z(updatedTraceback(k2, p2-3:p2-1)))];
                        vec5 = BB2 - BB1;
                        vec6 = BB1 - BB3;
                        vec7 = BB2 - BB4;
                        angle3 = acos(abs(dot(vec5, vec6)/(norm(vec5)*norm(vec6))));
                        angle3 = angle3/pi*180;
                        angle4 = acos(abs(dot(vec5, vec7)/(norm(vec5)*norm(vec7))));
                        angle4 = angle4/pi*180;
                        if angle3 < 40 && angle4 < 40
                            sig = true;
                        end
                    elseif p1 ~= NPR(k1) && p2~=1
                        BB1 = [x(updatedTraceback(k1, p1)), y(updatedTraceback(k1, p1)), z(updatedTraceback(k1, p1))];
                        BB2 = [x(updatedTraceback(k2, p2)), y(updatedTraceback(k2, p2)), z(updatedTraceback(k2, p2))];
                        BB3 = [x(updatedTraceback(k1, p1+1)), y(updatedTraceback(k1, p1+1)), z(updatedTraceback(k1, p1+1))];
                        BB4 = [x(updatedTraceback(k2, p2-1)), y(updatedTraceback(k2, p2-1)), z(updatedTraceback(k2, p2-1))];
                        vec5 = BB2 - BB1;
                        vec6 = BB1 - BB3;
                        vec7 = BB2 - BB4;
                        angle3 = acos(abs(dot(vec5, vec6)/(norm(vec5)*norm(vec6))));
                        angle3 = angle3/pi*180;
                        angle4 = acos(abs(dot(vec5, vec7)/(norm(vec5)*norm(vec7))));
                        angle4 = angle4/pi*180;
                        if angle3 < 40 && angle4 < 40
                            sig = true;
                        end
                    end
                        
                    if (angle2 < 5 || d2 < 3.5 || vecnorm(vec1, 2) < 2. || angle < 25 || sig) && checkOneRow(row, x, y, z)
%                         if p1 > 3
%                             updatedTraceback(end+1, 1:p1-1) = updatedTraceback(k1, 1:p1-1);
%                         end
%                         if NPR(k2)-p2 > 2
%                             updatedTraceback(end+1, 1:NPR(k2)-p2) = updatedTraceback(k2, p2+1:NPR(k2));
%                         end
                            
%                         updatedTraceback(k1, 1:NPR(k1)-p1+1) = updatedTraceback(k1, p1:NPR(k1));
%                         updatedTraceback(k1, NPR(k1)-p1+2:end) = 0;
%                         updatedTraceback(k1, NPR(k1)-p1+2:NPR(k1)-p1+p2+1) = updatedTraceback(k2, 1:p2);
                        updatedTraceback(k1, NPR(k1)+1:NPR(k1)+NPR(k2)) = updatedTraceback(k2,1:NPR(k2));
                        % update 'label' vector
                        updatedTraceback = sortByRow(updatedTraceback, k1, dist2Ant);
                        updatedTraceback(k2,:)=[];
                            
    %                     updatedTraceback = connect1Line(updatedTraceback, k1, k2, NPR(k1), NPR(k2), dist2Ant, dist2All);
                        changed = true;
                    end
                end
            end
        end
        if changed
            break
        end
    end
end
end