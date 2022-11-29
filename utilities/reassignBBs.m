% -------------------------------------------------------------------------
% [Ben] 2/2/18 (added by Ben)
% Performs the greedy reassignment of BBs to most appropriate ciliary rows.
% Works on an updatedTraceback matrix that already only contains legal
% rows (the output of getIniLabel). Returns the F-score of the final
% ciliary row configuration arrived at, as well as the corresponding
% traceback matrix.
% -------------------------------------------------------------------------


function [f_final, initialMatrix] = reassignBBs(updatedTraceback, x, y, z, antPole, postPole, dist2Ant, dist2Post, withLabel,scale_xyz)

coeffs = getWeights(updatedTraceback, antPole, postPole, x, y, z, dist2Ant, dist2Post,scale_xyz);
coeffs = [0.3, 0.3, 1, 0, 1, 0].*coeffs;
% coeffs = [0.33, 0.11, 0.11, 0.11, 0, 0.33].*coeffs;
% coeffs = [0.2, 0.2, 0.2, 0.2, 0.2, 0.2].*coeffs;
%value of objective function
[f_old, ~] = constraints(updatedTraceback, antPole, postPole, dist2Ant, dist2Post, x, y, z, coeffs,scale_xyz);

initialMatrix = updatedTraceback;
[numRow, ~] = size(initialMatrix);
numBB_withLabel = length(withLabel);

% fprintf('numRow %d numBB_withLabel %d\n',numRow, numBB_withLabel);

% pg. 20 of Jingyi's paper. Mtx is a matrix with n rows and m columns,
% where n is the no. of valid ciliary rows in the cell, and m is the no. of
% BBs that have been assigned to a ciliary row. Mtx(a, b) saves the
% difference between the F-score of the current BB assignment and the
% assignment if we assigned BB 'b' to row 'a'.
minMtxVal = -Inf;

% plot constraints
figure(1)
hold on
c1 = animatedline('marker', 'o');
c2 = animatedline('marker', '+');
c3 = animatedline('marker', '*');
c4 = animatedline('marker', '.');
c5 = animatedline('marker', 'x');
c6 = animatedline('marker', 'd');
ct = animatedline('color', 'r');
iter = 0;
while minMtxVal < 0
    iter = iter + 1;
    Mtx = zeros(numRow, numBB_withLabel);
    for row = 1:numRow
        parfor bb = 1:numBB_withLabel
            %newMtx = moveBB(initialMatrix, withLabel(bb), row, dist2Ant,scale_xyz);
            newMtx = moveBB(initialMatrix, withLabel(bb), row, dist2Ant);
            [f, ~] = constraints(newMtx, antPole, postPole, dist2Ant, dist2Post, ...
                            x, y, z, coeffs,scale_xyz);
            Mtx(row, bb) = f - f_old;
            
        end
    end
    
    % do not suppress output so we can observe the convergence 
    minMtxVal = min(Mtx,[],'all')
    f_old = f_old + minMtxVal;
    
    % if there are several entries that all have the same min value, then
    % just pick the first pair of indices
    [minValRowIdx, minValBBIdx] = find(Mtx == minMtxVal, 1);
    initialMatrix = moveBB(initialMatrix, withLabel(minValBBIdx), minValRowIdx, dist2Ant);
    
    % to show the constraints
    [f, fs] = constraints(initialMatrix, antPole, postPole, dist2Ant, dist2Post, ...
                            x, y, z, coeffs,scale_xyz);
                        
    % draw CR in each iteration
    fig = figure(3);
    clf(fig,'reset')
    % visualize_bbs(x, y, z, initialMatrix, antPole, postPole, scale_xyz);              

    addpoints(c1,iter,fs(1))
    addpoints(c2,iter,fs(2))
    addpoints(c3,iter,fs(3))
    addpoints(c4,iter,fs(4))
    addpoints(c5,iter,fs(5))
    addpoints(c6,iter,fs(6))
    addpoints(ct,iter,f)
    legend([c1 c2 c3 c4 c5 c6, ct],{'c1','c2','c3','c4','c5','c6', 'Total'})
    drawnow
 
end

[~,y1] = getpoints(c1);
[~,y2] = getpoints(c2);
[~,y3] = getpoints(c3);
[~,y4] = getpoints(c4);
[~,y5] = getpoints(c5);
[~,y6] = getpoints(c6);
[~,yt] = getpoints(ct);
save('1005_161101.mat','y1','y2','y3','y4','y5','y6','yt', 'coeffs');
f_final = f_old;
end

