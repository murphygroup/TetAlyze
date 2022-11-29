clc
clear all
close all
main_path = '2- Cell cycle generative model'; 
paths = ["N=1"; "N=2"; "N=3";];

data = [];
skipped = 0;
ratio = [];
nBB = [];
D = dir(main_path); 
for i=1:length(paths)
    full_path = fullfile(main_path, paths(i));
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.
    for k = 3:length(D) % avoid using the first ones
        currD = D(k).name;
        tempPath = fullfile(full_path, currD);
        if isfolder(tempPath)
            fList = dir(tempPath); % Get the file list in the subdirectory
            for j = 3:length(fList)
                currD = fList(j).name;
                if contains(currD, '.fig') && ~contains(currD, '._')
                    [p, waist_, n] = measurePerimeter(fullfile(tempPath, currD));
                       nBB(end + 1) = n;
                       ratio(end + 1) = waist_/p;
                    break
                end
            end
        end
    end
end
s = scatter(nBB, ratio, 'MarkerEdgeColor', '#D95319', 'MarkerFaceColor', '#D95319');
s.AlphaData = 0.5;
s.MarkerFaceAlpha = 0.5;
s.MarkerEdgeAlpha = 0.5;
xlabel('Number of BBs', 'FontSize', 14);
ylabel('middle region perimeter over cell perimeter', 'FontSize', 14);
xlim([350 750]);
xticks(350:50:750);
ax = gca; 
ax.FontSize = 14

function [perim waist n] = measurePerimeter(figname)
    if (exist(figname,'file'))
        h1=openfig(figname, 'invisible');
    end
    ax = gca; 
    h = findobj(gca,'Type','line'); 

    x = [];
    y = [];
    z = [];
    
    rowIdx = [];
    for i=1:length(h)   
        x = [x h(i).XData];
        y = [y h(i).YData];
        z = [z h(i).ZData];
        for k=1:length(h(i).ZData)
            rowIdx(end + 1) = i;
        end
    end
    
    close(h1)
    minZ = min(z);
    maxZ = max(z);
    
    xMid = [];
    yMid = [];
    rowIdx2 = [];
    for i=1:length(x)
        if z(i) > minZ + 0.45 * (maxZ - minZ) & z(i) < minZ + 0.55 * (maxZ - minZ)
            xMid(end + 1) = x(i);
            yMid(end + 1) = y(i);
            rowIdx2(end + 1) = rowIdx(i);
        end
    end
    data = vertcat(x,y)';
    [coeff,score,latent] = pca(data);
    k = convhull(score(:,1:2));
%     figure()
    perim = 0;
    convRowIdx = rowIdx(k);
    for i=1:length(k)-1
        perim = perim + norm(score(k(i),1:2) - score(k(i+1),1:2));
%         plot(score(k(i:i+1),1),score(k(i:i+1),2));
%         hold on
    end
    
    data = [xMid' yMid'];
    n = length(x);
    score2 = data*coeff;
    k = convhull(score2(:,1:2));
	convRowIdx2 = rowIdx2(k);
    waist = 0;
    for i=1:length(k)-1
        waist = waist + norm(score2(k(i),1:2) - score2(k(i+1),1:2));
%         plot(score2(k(i:i+1),1),score2(k(i:i+1),2));
%         hold on
    end
    r = waist/perim;
    if waist/perim < 0.93
        figure();
        k = convhull(score(:,1:2));
        for i=1:length(k)-1
            if i == length(k)-1
                plot(score(k(i:i+1),1),score(k(i:i+1),2), 'o-', 'LineWidth', 2, 'Color', 'black');
            else
                plot(score(k(i:i+1),1),score(k(i:i+1),2), 'LineWidth', 2, 'Color', 'black');
            end
            hold on
        end
        k = convhull(score2(:,1:2));
        for i=1:length(k)-1
            plot(score2(k(i:i+1),1),score2(k(i:i+1),2), 'LineWidth', 2, 'Color', 'red');
        end
        axis equal
        axis off
        title(sprintf('BB#=%d', n))
    end
end
