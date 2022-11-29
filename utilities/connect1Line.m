function updatedTraceback = connect1Line(updatedTraceback, k1, k2, rowNum1, rowNum2, dist2Ant, dist2All)

updatedTraceback(k1, rowNum1+1:rowNum1+rowNum2) = updatedTraceback(k2,1:rowNum2);
% update 'label' vector
updatedTraceback = sortByRow(updatedTraceback, k1, dist2Ant);
updatedTraceback(k2,:)=[];
% row1 = updatedTraceback(k1, 1:rowNum1);
% row2 = updatedTraceback(k2, 1:rowNum2);
% localDist = dist2All(row1, row2);
% minDist = min(localDist(:));
% [p1, p2] = find(localDist == minDist);
% updatedTraceback(k1, 1:rowNum1-p1+1) = updatedTraceback(k1, p1:rowNum1);
% updatedTraceback(k1, rowNum1-p1+1:rowNum1-p1+p2) = updatedTraceback(k2, 1:p2);
% % update 'label' vector
% updatedTraceback = sortByRow(updatedTraceback, k1, dist2Ant);
% updatedTraceback(k2,:)=[];
end
                