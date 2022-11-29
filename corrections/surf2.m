function surf(BW, cort_x, cort_y, cort_z, antIdx, postIdx, midIdx, resultPath)


P = alphaShape(cort_x,cort_y,cort_z);
P.Alpha=20;
% figure(10)
% plot(P)
antP = alphaShape(cort_x(antIdx),cort_y(antIdx),cort_z(antIdx));
% figure(11)
% plot(antP)
antP.Alpha = 20;
antArea = surfaceArea(antP);

entireArea = surfaceArea(P);

notAntIdx = setdiff(1:length(BW), antIdx);
notAntP = alphaShape(cort_x(notAntIdx),cort_y(notAntIdx),cort_z(notAntIdx));
notAntP.Alpha = 20;
% figure(12)
% plot(notAntP)
notAntArea = surfaceArea(notAntP);
tempArea1 = (notAntArea + antArea - entireArea)/2;


postP = alphaShape(cort_x(postIdx),cort_y(postIdx),cort_z(postIdx));
postP.Alpha = 20;
% figure(13)
% plot(postP)
postArea = surfaceArea(postP);

notPostIdx = setdiff(1:length(BW), postIdx);
notPostP = alphaShape(cort_x(notPostIdx),cort_y(notPostIdx),cort_z(notPostIdx));
notPostP.Alpha = 20;
% figure(14)
% plot(notPostP)

notPostArea = surfaceArea(notPostP);
tempArea2 = (notPostArea + postArea - entireArea)/2;

antArea = antArea - tempArea1;
postArea = postArea - tempArea2;

midP = alphaShape(cort_x(midIdx),cort_y(midIdx),cort_z(midIdx));
midP.Alpha = 20;
midArea1 = surfaceArea(midP) - tempArea1 - tempArea2;

fid = fopen(join([resultPath, '/', 'Summary.txt']),'a');
fprintf(fid, 'Anterior Area (um^2): %f\n', antArea);
fprintf(fid, '# Anterior BBs: %d\n', length(antIdx));

fprintf(fid, 'Medial Area (um^2): %f\n', midArea1);
fprintf(fid, '# Medial BBs: %d\n', length(midIdx));

fprintf(fid, 'Posterior Area (um^2): %f\n', postArea);
fprintf(fid, '# Posterior BBs: %d\n', length(postIdx));

fprintf('Anterior Area (um^2): %f\n', antArea);
fprintf('# Anterior BBs: %d\n', length(antIdx));

fprintf('Medial Area (um^2): %f\n', midArea1);
fprintf('# Medial BBs: %d\n', length(midIdx));

fprintf('Posterior Area (um^2): %f\n', postArea);
fprintf('# Posterior BBs: %d\n', length(postIdx));

end