function [startPt, endPt] = findClosest(candidates, distances, heads, tails)
[d1, ~] = size(distances);
startPt = [];
endPt = [];
closest = zeros(d1, 1);

for i=1:d1
    dists = distances(i,:); % length: len(pts_end)
    [~, idx] = sort(dists, 'ascend'); % gives ordering to rearrange in ascending order
    if ismember(candidates(i), heads)
        k = find(heads == candidates(i));
        if candidates(idx(2)) == tails(k)
            closest(i) = idx(3);
        else
            closest(i) = idx(2);
        end   
    elseif ismember(candidates(i), tails)
        k = find(tails == candidates(i));
        if candidates(idx(2)) == heads(k)
            closest(i) = idx(3);
        else
            closest(i) = idx(2);
        end
    else
        closest(i) = idx(2);
    end
end
for i=1:d1
    j = closest(i);
    if i < j
       k = closest(j);
       if k == i
           startPt(end+1) = candidates(i);
           endPt(end+1) = candidates(j);
       end
    end
end

end