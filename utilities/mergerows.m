% -------------------------------------------------------------------------
% [Chongming] 8/29/20
% merge close rows
% -------------------------------------------------------------------------
function [merged_matrix , mergeMatrix] = mergerows(x, y, z, mtx, antPole, postPole, dist2Ant)

mtx = sort_d2ant(mtx, dist2Ant);

row_counts = sum(mtx ~= 0, 2);
num_rows = length(row_counts);
distThres = 30;
angleThres1 = 30/180*pi;
angleThres2 = 15/180*pi;
ant2post = postPole - antPole;
mergeMatrix = zeros(num_rows);

% merge two rows if their distance is less than distThres and satisfy angle
% constraint.
for i = 1:num_rows
    istart = [x(mtx(i, 1)) y(mtx(i, 1)) z(mtx(i, 1))];
    iend = [x(mtx(i, row_counts(i))) y(mtx(i, row_counts(i))) z(mtx(i, row_counts(i)))];
    veci = istart - iend;
    for j = i+1:num_rows
        jstart = [x(mtx(j, 1)) y(mtx(j, 1)) z(mtx(j, 1))];
        jend = [x(mtx(j, row_counts(j))) y(mtx(j, row_counts(j))) z(mtx(j, row_counts(j)))];
        vecj = jstart - jend;
        vec1 = iend- jstart;
        vec2 = jend - istart;
        angle1 = acos(dot(ant2post, vec1)/(norm(ant2post)*norm(vec1)));
        angle2 = acos(dot(ant2post, vec2)/(norm(ant2post)*norm(vec2)));
        anglei1 = acos(dot(veci, vec1)/(norm(veci)*norm(vec1)));
        anglei2 = acos(dot(veci, vec2)/(norm(veci)*norm(vec2)));
        anglej1 = acos(dot(vecj, vec1)/(norm(vecj)*norm(vec1)));
        anglej2 = acos(dot(vecj, vec2)/(norm(vecj)*norm(vec2)));
        %anglei1 = min(anglei1, pi-anglei1);
        %anglei2 = min(anglei2, pi-anglei2);
        %anglej1 = min(anglej1, pi-anglej1);
        %anglej2 = min(anglej2, pi-anglej2);
        %angle1 = min(angle1, pi-angle1);
        %angle2 = min(angle2, pi-angle2);
        %ij_dist = min(norm(istart-jend), norm(iend-jstart));
        if (max(anglei1, anglej1) < angleThres2 && norm(vec1) < distThres && angle1 < angleThres1) || ...
           (max(anglei2, anglej2) < angleThres2 && norm(vec2) < distThres && angle2 < angleThres1)
            mergeMatrix(i, j) = 1;
        end
    end
end

% merge rows using union-find data structure
dj = DJSet(num_rows);
for i = 1:num_rows
    for j = i+1:num_rows
        if mergeMatrix(i, j) == 1
            dj.union(i, j)
        end
    end
end

root_rows_idx = zeros(num_rows, 1);
root_rows_counts = zeros(num_rows, 1);
for i = 1:num_rows
    root_rows_idx(dj.find(i)) = 1;
    root_rows_counts(dj.find(i)) = root_rows_counts(dj.find(i)) + row_counts(i);
end
max_root_rows_counts = max(root_rows_counts);
merged_matrix = zeros(num_rows, max_root_rows_counts);
merged_matrix_idx = zeros(num_rows, 1);
for i = 1:num_rows
    root = dj.find(i);
    root_idx = merged_matrix_idx(root)+1;
    merged_matrix(root, root_idx:root_idx+row_counts(i)-1) = mtx(i, 1:row_counts(i));
    merged_matrix_idx(root) = root_idx + row_counts(i)-1;
end

merged_matrix = merged_matrix(root_rows_idx>0, :);
merged_matrix = sort_d2ant(merged_matrix, dist2Ant);

end