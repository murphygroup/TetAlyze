function BBstat(I, BW, cort_x, cort_y, cort_z, resultPath, Iori)


pxIntensities = zeros(length(BW), 1);
pxIntensitiesOri = zeros(length(BW), 1);
a = length(BW);
zz = cast(BW, 'uint32');
antIdx = [];
postIdx = [];
midIdx = [];
a = max(cort_z) ;
b = min(cort_z);
for i=1:length(BW)
    pxIntensities(i) = I(zz(i, 2), zz(i, 1), zz(i, 3));
    pxIntensitiesOri(i) = Iori(zz(i, 2), zz(i, 1), zz(i, 3));
    if cort_z(i) >= a - 0.25*(a-b)
        antIdx(end+1) = i;
    elseif cort_z(i) <= b + 0.25*(a-b)
        postIdx(end+1) = i;
    else
        midIdx(end+1) = i;
    end
end
writematrix(pxIntensitiesOri, join([resultPath, '/', 'Raw BB intensity (all).csv']));
writematrix(pxIntensitiesOri(antIdx), join([resultPath, '/', 'Raw BB intensity (anterior).csv']));
writematrix(pxIntensitiesOri(midIdx), join([resultPath, '/', 'Raw BB intensity (medial).csv']));
writematrix(pxIntensitiesOri(postIdx), join([resultPath, '/', 'Raw BB intensity (posterior).csv']));
f = figure(5);
f1 = subplot(2,2,1);

antAvg = mean(pxIntensities(antIdx));
pxIntensities = pxIntensities / antAvg;
maxValue = max(pxIntensities)+0.5;
h1 = histogram(pxIntensities, 20, 'Normalization', 'probability');
title('All BBs intensities (intensities are normalized to 0-1)')
xlim([0 maxValue]); 
xlabel('Intensity')
ylabel('BBs frequency')

f2 = subplot(2,2,2);
h2 = histogram(pxIntensities(antIdx), 20, 'Normalization', 'probability');
title('Anterior BBs intensities')
xlim([0 maxValue]); 
xlabel('Intensity')
ylabel('BBs frequency')


f3 = subplot(2,2,3);
h3 = histogram(pxIntensities(midIdx), 20, 'Normalization', 'probability');
title('Medial BBs intensities')
xlim([0 maxValue]); 
xlabel('Intensity')
ylabel('BBs frequency')

f4 = subplot(2,2,4);
h4 = histogram(pxIntensities(postIdx), 20, 'Normalization', 'probability');
title('Posterior BBs intensities')
xlim([0 maxValue]); 
xlabel('Intensity')
ylabel('BBs frequency')
saveas(f, join([resultPath, '/', 'BB Intensity Histogram.png']));

surf_local(BW, cort_x, cort_y, cort_z, antIdx, postIdx, midIdx, resultPath, 'Summary.txt');

yMax = max([max(h1.Values), max(h2.Values), max(h3.Values), max(h4.Values)]);
set(f1, 'YLim', [0 yMax]);
set(f2, 'YLim', [0 yMax]);
set(f3, 'YLim', [0 yMax]);
set(f4, 'YLim', [0 yMax]);
% dist2OA = vecnorm([cort_x cort_y cort_z] - lowestOA(1, :),2,2);
% [~, order2OA] = sort(dist2OA,'ascend');


end