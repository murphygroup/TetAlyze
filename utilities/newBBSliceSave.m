function newBBSliceSave(BW, I, resultPath, imageID)

imgStack = [];
startFrame = 1000000;
endFrame = -1;
z = cast(BW(:, 3), 'uint32');
[height, width, zLength] = size(I);

for i=1:zLength
    idx = find(z == i);
    Xs = BW(idx, 1);
    Ys = BW(idx, 2);
    if ~isempty(idx)
        if i < startFrame
            startFrame = i;
        end
        if i > endFrame
            endFrame = i;
        end
        img = I(:, :, i)*10;
        img = cat(3, img, img, img);
        imgStack = cat(4, imgStack, img);
        for j =1:length(idx)
            img = insertMarker(img,[Xs(j) Ys(j)],'circle', 'size',4);
        end

        imgStack = cat(4, imgStack, img);
    end
end
if ~isempty(imgStack)
    imgStack = permute(imgStack, [1 2 4 3]);
    % a = size(imgStack)
    figure(10);
    sliceViewer(imgStack);

    C = strsplit(resultPath, '\');
%     save(char(fullfile('E:\BB Project\big-1\new BB slices', strcat(C(end), ' new BB Identification.mat'))), 'imgStack');
end
end