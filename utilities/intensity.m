clc
close all
clear all


main_path = 'E:\BB Project\generative model images\2- Cell cycle generative model';   %  Path to the main results folder
thetaAll = [];
thetaAllMid = [];
phiAll = [];
colors = [];

data_ = table2array(readtable('WT-new.csv'));

paths = ["N=1"; "N=2"; "N=3";];
% paths = ["B1868 t=8"];
values = zeros(338, 4);
j1 = 1;
j2 = 1;
j3 = 1;
j4 = 1;
j = 1;
r = 2;
N = 20;
count = zeros(N, 1);
count1 = zeros(N, 4);
count2 = zeros(N, 4);
count3 = zeros(N, 4);
count4 = zeros(N, 4);
midCount = [];
spacing = [];

for i=1:length(paths)
    full_path = fullfile(main_path, paths(i));
%         full_path = 'E:\BB Project\4- B1868_ _Poc1-mCh cell division';
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.

    for k = 3:length(D) % avoid using the first ones
        if D(k).isdir
            currD = D(k).name; % Get the current subdirectory name

            tempPath = fullfile(full_path, currD);

            data = csvread(fullfile(tempPath, 'intensityAnalysis.csv'));
            alignment = csvread(fullfile(tempPath, 'Alignment.csv'));
            data2 = csvread(fullfile(tempPath, 'data4PCA.csv'));
            
%             if length(data) > 400 
%                 continue
%             end
            if size(alignment, 1) ~= N
                continue
            end
% 
% %             
            thetaAll = [thetaAll data(:, 1)'];
            midIdx = [];
            phiAll = [phiAll data(:, 2)'];
            colors = [colors (data(:, 4)>r)'];
%             
%             
%             idx_ =  find(data(:, 4)>r);
            countTemp = zeros(size(alignment, 1), 1);
            for ii=1:length(data)
                [row ~]=find(alignment == ii);
                if -0.5 < cos(data(ii, 2)) && cos(data(ii, 2)) < 0.5
                    countTemp(row) = countTemp(row) + 1;
                end
%                 [row ~]=find(alignment == idx_(ii));
%                 if cos(data(idx_(ii), 2)) > 0.5
%                     countTemp(row, 1) = countTemp(row, 1) + 1;
%                 elseif 0 < cos(data(idx_(ii), 2)) & cos(data(idx_(ii), 2)) < 0.5
%                     countTemp(row, 2) = countTemp(row, 2) + 1;
%                 elseif -0.5 < cos(data(idx_(ii), 2)) & cos(data(idx_(ii), 2)) < 0
%                     countTemp(row, 3) = countTemp(row, 3) + 1;
%                 else
%                     countTemp(row, 4) = countTemp(row, 4) + 1;
%                 end
            end
%             
%             
%             midCount(end + 1) = sum(countTemp(:, 2)) + sum(countTemp(:, 3));
%             spacing(end + 1) = data2(15);
%             countTemp = countTemp / length(idx_);
            count = count + countTemp;
%             if length(data) <= 400
%                 count1 = count1 + countTemp;
%                 j1 = j1 + 1;
%             elseif length(data) > 400 & length(data) <= 500
%                 count2 = count2 + countTemp;
%                 j2 = j2 + 1;
%             elseif length(data) > 500 & length(data) <= 600
%                 count3 = count3 + countTemp;
%                 j3 = j3 + 1;
%             else
%                 count4 = count4 + countTemp;
%                 j4 = j4 + 1;
%             end
%             

                
            
            
        j = j + 1;    
        end
                
    end

end

count = count / j;
count1 = count1 / j1;
count2 = count2 / j2;
count3 = count3 / j3;
count4 = count4 / j4;
idx = find(colors);



antIdx = find(cos(phiAll)>0.5);
nAnt = length(nonzeros(ismember(idx, antIdx)))/length(antIdx);

postIdx = find(cos(phiAll)<-0.5); 
nPost = length(nonzeros(ismember(idx, postIdx)))/length(postIdx);
% 
medIdx1 = find(-0.5<cos(phiAll) & cos(phiAll)<0); 
nMed1 = length(nonzeros(ismember(idx, medIdx1)))/length(medIdx1);

medIdx2 = find(0<cos(phiAll) & cos(phiAll)<0.5); 
nMed2 = length(nonzeros(ismember(idx, medIdx2)))/length(medIdx2);


values(j, 1) = nAnt;
values(j, 2) = nMed2;
values(j, 3) = nMed1;
values(j, 4) = nPost;

nMed = (length(idx) - length(nonzeros(ismember(idx, antIdx))) - length(nonzeros(ismember(idx, postIdx))))/(length(thetaAll')-length(antIdx)-length(postIdx));



figure();
t = tiledlayout(5,1,'TileSpacing','Compact','Padding','Compact');
nexttile;
xvalues = {'row 1', 'row 2', 'row 3', 'row 4', 'row 5', 'row 6', 'row 7', 'row 8', 'row 9', 'row 10',...
    'row 11', 'row 12', 'row 13', 'row 14', 'row 15', 'row 16', 'row 17', 'row 18', 'row 19'};
yvalues = {'anterior','anterior medial','posterior medial','posterior'};
heatmap(xvalues, yvalues, count');
title('N = 19, combined');
nexttile;
xvalues = {'row 1', 'row 2', 'row 3', 'row 4', 'row 5', 'row 6', 'row 7', 'row 8', 'row 9', 'row 10',...
    'row 11', 'row 12', 'row 13', 'row 14', 'row 15', 'row 16', 'row 17', 'row 18', 'row 19'};
yvalues = {'anterior','anterior medial','posterior medial','posterior'};
heatmap(xvalues, yvalues, count1');
title('BB# <= 400');
nexttile;
xvalues = {'row 1', 'row 2', 'row 3', 'row 4', 'row 5', 'row 6', 'row 7', 'row 8', 'row 9', 'row 10',...
    'row 11', 'row 12', 'row 13', 'row 14', 'row 15', 'row 16', 'row 17', 'row 18', 'row 19'};
yvalues = {'anterior','anterior medial','posterior medial','posterior'};
heatmap(xvalues, yvalues, count2');
title('400 < BB# <= 500');
nexttile;
xvalues = {'row 1', 'row 2', 'row 3', 'row 4', 'row 5', 'row 6', 'row 7', 'row 8', 'row 9', 'row 10',...
    'row 11', 'row 12', 'row 13', 'row 14', 'row 15', 'row 16', 'row 17', 'row 18', 'row 19'};
yvalues = {'anterior','anterior medial','posterior medial','posterior'};
heatmap(xvalues, yvalues, count3');
title('500 < BB# <= 600');
nexttile;
xvalues = {'row 1', 'row 2', 'row 3', 'row 4', 'row 5', 'row 6', 'row 7', 'row 8', 'row 9', 'row 10',...
    'row 11', 'row 12', 'row 13', 'row 14', 'row 15', 'row 16', 'row 17', 'row 18', 'row 19'};
yvalues = {'anterior','anterior medial','posterior medial','posterior'};
heatmap(xvalues, yvalues, count4');
title('BB# > 600');
% scatter3(cos(thetaAll(idx)).*sin(phiAll(idx)), sin(thetaAll(idx)).*sin(phiAll(idx)), cos(phiAll(idx)), 5, colors(idx), 'red');
[x, y, z] = sphere(50);
% [xNew,yNew,zNew] = cart2sph(x(:), y(:), z(:));
y_b = mvksdensity([(cos(thetaAll(idx)).*sin(phiAll(idx)))' (sin(thetaAll(idx)).*sin(phiAll(idx)))' cos(phiAll(idx))'], [x(:), y(:), z(:)], 'Bandwidth', 0.25, 'Kernel', 'normal');
y_b2 = mvksdensity([(cos(thetaAll).*sin(phiAll))' (sin(thetaAll).*sin(phiAll))' cos(phiAll)'], [x(:), y(:), z(:)], 'Bandwidth', 0.25, 'Kernel', 'normal');
figure();
s = surf(x, y, z, reshape(y_b, [51, 51]));
colormap(jet);
s.EdgeColor = 'none';
title('spatial desnsity map of immature cortical BBs')
axis off
figure();
s = surf(x, y, z, reshape(y_b2, [51, 51]));
colormap(jet);
s.EdgeColor = 'none';
title('spatial desnsity map of cortical BBs')
axis off
figure();
scatter(spacing, midCount, 'o');
r = corrcoef(spacing, midCount);
text(1.6, 75, sprintf('Pearson correlation = %0.3f', r(2, 1)), 'FontSize', 15);
xlabel('average BB spacing');
ylabel('new BB #');

% colorbar;
% plot3(x(1, :), x(:, 1), z(:, 1));
% R1 = corrcoef(data_(:, 1), values(:, 1));
% R2 = corrcoef(data_(:, 2), values(:, 2));
% R3 = corrcoef(data_(:, 3), values(:, 3));
% R4 = corrcoef(data_(:, 4), values(:, 4));

% figure();
% plot(data_(:, 1), values(:, 1), '.', data_(:, 1),values(:, 2), '.', data_(:, 1),values(:, 3), '.', data_(:, 1),values(:, 4), '.');

% legend('anterior','anterior medial','posterior medial','posterior');
% xlabel('ratio');
% ylabel('fraction');
% title(paths);