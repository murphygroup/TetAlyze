% -------------------------------------------------------------------------
% [Ben] 12/07/17
% Returns BBs that are part of a valid anterior-posterior pair, together
% with their corresponding partner (i.e. the other half making up the pair)
% ------------------------------------------------------------------------- 


function [startPt, endPt] = findLink(x, y, z, distances, antPole, postPole, dist2Ant, dist2Post)
n = length(x);
pts_start = 1:n;   % start: posterior
pts_end = 1:n;    % end: anterior

% add distance from point i to plane defined by (j,antPole,postPole) to
% pairwise distance between point i and j
% min dist between partner and plane(BB maxima, antPole, postPole)

% So dist_forward seems to be a matrix where entry (i,j) stores the sum of
% the distance between BBs i and j (indexed based on the order of
% getBBIdx's output) and the distance between j and the plane defined by
% the 2 poles (findPoles.m) and i.

% This metric calculated in both the forward and reverse directions and BBs
% are assigned a partner only if the forward and reverse directions
% identify the same reciprocal BB pair.


dist_backward = distances;
dist_forward = distances;
% You find the anterior partner of i by finding the j such that
% dist_forward(i,j) is smallest and dist2Ant(j) < dist2Ant(i).
% You find the posterior parter of i by finding the j such that
% dist_backward(i,j) is smallest and dist2Ant(j) > dist2Ant(i).

startPt = [];
endPt = [];
prevLen = n;

startPt_forward = zeros(n ,1);
endPt_forward = zeros(n, 1);
% finds anterior partners using forward direction results
for i = 1:n
    distances = dist_forward(i,:); % length: len(pts_end)
    [~, idx] = sort(distances, 'ascend'); % gives ordering to rearrange in ascending order
    for j = 1:n
        if (idx(j) ~= i) && ...
                (dist2Ant(idx(j)) < dist2Ant(i)) % partner is anterior neighbor
            startPt_forward(i) = i; % why not just startPt_forward(i) = i?
            endPt_forward(i) = idx(j); 
            break
        end
    end
end

startPt_backward=zeros(n, 1);
endPt_backward=zeros(n, 1);
% finds posterior partners using backward direction results
for i = 1:n
    distances = dist_backward(i, :); % length: len(pts_start)
    [~, idx] = sort(distances, 'ascend');
    for j = 1:n
        if (idx(j) ~= i) && ...
                (dist2Ant(idx(j)) > dist2Ant(i)) % partner is posterior neighbor
            startPt_backward(i) = i;
            endPt_backward(i) = idx(j);
            break
        end
    end
end

% So now endPt_forward contains all the anterior neighbors if they exist,
% and endPt_backward contains all the posterior neighbors if they exist.
% Note that startPt_forward and startPt_backward will have a 0 at entry i
% if BB_i does not have an anterior partner and does not have a posterior
% partner respectively.

% At the end of the connection process, the only unconnected BBs that
% remain are those that will never yield reciprocal connections. These BBs
% are connected to the closest unconnected BB only if the connection is
% less than 5um. 
if ~isempty(startPt_forward)
    for i = 1:length(startPt_forward)
        % BB_i has anterior partner
        if startPt_forward(i) ~= 0 && endPt_forward(i) ~= 0
            % find anterior partner of BB_i
            idx = find(startPt_backward == endPt_forward(i));
            % check that posterior partner of bb_i's anterior partner is in
            % fact BB_i
            if endPt_backward(idx) ~= 0 && endPt_backward(idx) == startPt_forward(i)
                % p1 is BB_i and p2 is p1's anterior partner
                p1 = [x(startPt_forward(i)), y(startPt_forward(i)), z(startPt_forward(i))];
                p2 = [x(endPt_forward(i)), y(endPt_forward(i)), z(endPt_forward(i))];
                % distance between BB and its partner < 5um
                
                % -------------------------------------------------------------------------
                % [Chongming] 08/29/20
                % Add angle constraints to avoid zig-zag or protrusion
                % -------------------------------------------------------------------------
                ant2post = postPole - antPole;
                p22p1 = p1 - p2;
                angleWithPoles = fourptsangle(p1, p2, antPole, postPole);
                
                if angleWithPoles < 35/180*pi 
                    % startPt and endPt together record pairs of BBs and
                    % their anterior partners.
                    startPt = [startPt, startPt_forward(i)];
                    endPt = [endPt, endPt_forward(i)];
                end
                
            end
        end
    end
end

end






