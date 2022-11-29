clc
clear all
N = 19;
N2= 20;
data_ = table2array(readtable('WT-new.csv'));
main_path = '2- Cell cycle generative model'; 
paths = ["N=1"; "N=2"; "N=3";];
midCount = [];
spacing = zeros(300, N, 4);
nCell = 1;
nCell2 = 1;
skipped = 0;
r = 2;
counts = zeros(300, N, 4);
centers = zeros(N, 1);
centers2 = zeros(N2, 1);

spacing = zeros(N, 1);
spacing2 = zeros(N2, 1);


D = dir(main_path); 
for j=1:length(paths)
    full_path = fullfile(main_path, paths(j));
%         full_path = 'E:\BB Project\4- B1868_ _Poc1-mCh cell division';
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.

    for k = 3:length(D) % avoid using the first ones
        if D(k).isdir
            currD = D(k).name; % Get the current subdirectory name

            tempPath = fullfile(full_path, currD);

            
            data = csvread(fullfile(tempPath, 'intensityAnalysis.csv'));
            
            alignment = csvread(fullfile(tempPath, 'Alignment.csv'));
            spacing_data = csvread(fullfile(tempPath, 'BB Row Spacing Raw Data.csv'));
            if sum(isnan(spacing_data(:, 3))) > 0
                continue
            end
            if size(alignment, 1) == N
                nRow = size(alignment, 1);
                center_temp = zeros(N, 1);
                theta = data(:, 1);
                phi = data(:, 2);
%                 figure(8);

                for ii=1:nRow
                    idx_temp = [];
                    for jj=1:length(nonzeros(alignment(ii, :)))
                        if true
                            idx_temp(end+1) = alignment(ii, jj);
                        end
                    end
                    center_temp(ii) = std(data(idx_temp, 1));
                    if center_temp(ii) < 0
                        center_temp(ii) = 2*pi + center_temp(ii);
                    end
                    if center_temp(ii) > 2*pi
                        center_temp(ii) = center_temp(ii) - 2*pi;
                    end
%                     plot3(sin(phi(idx_temp)).*cos(theta(idx_temp)), sin(phi(idx_temp)).*sin(theta(idx_temp)), cos(phi(idx_temp)), 'o-', 'MarkerSize', 6);
%                     axis equal
%                     hold on
                end

                centers = centers + center_temp;
                nCell = nCell + 1;
                spacing = spacing + spacing_data(:, 3);
            elseif size(alignment, 1) == N2
                nRow = size(alignment, 1);
                center_temp = zeros(N2, 1);
                theta = data(:, 1);
                phi = data(:, 2);
    %             figure(8);

                for ii=1:nRow
                    idx_temp = [];
                    for jj=1:length(nonzeros(alignment(ii, :)))

                        if true
                            idx_temp(end+1) = alignment(ii, jj);
                        end
                    end
                    center_temp(ii) = std(data(idx_temp, 1));
                    if center_temp(ii) < 0
                        center_temp(ii) = 2*pi + center_temp(ii);
                    end
                    if center_temp(ii) > 2*pi
                        center_temp(ii) = center_temp(ii) - 2*pi;
                    end
    %                 scatter(sin(phi(idx_temp)).*cos(theta(idx_temp)), sin(phi(idx_temp)).*sin(theta(idx_temp)), 20 );
    %                 hold on
                end

                centers2 = centers2 + center_temp;
                nCell2 = nCell2 + 1;
                spacing2 = spacing2 + spacing_data(:, 3);
            end
        end
    end
end

centers = centers / nCell;
centers2 = centers2 / nCell2;
spacing = spacing / nCell;
spacing2 = spacing2 / nCell2;
figure();
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');

nexttile
plot(centers, 'o-', 'LineWidth', 1.5);
hold on
plot(centers2, 'o-', 'LineWidth', 1.5);
xlabel('Row index', 'FontSize', 14);
legend('19 ciliary rows','20 ciliary rows', 'FontSize', 14);
ylabel('Std(\theta)', 'FontSize', 14);

title('Ciliary row twistedness', 'FontSize', 14);

nexttile
plot(spacing, 'o-', 'LineWidth', 1.5);
hold on
plot(spacing2, 'o-', 'LineWidth', 1.5);


xticks(1:2:20);
xlabel('Row index', 'FontSize', 14);
legend('19 ciliary rows','20 ciliary rows', 'FontSize', 14);
ylabel('Row spacing (\mum))', 'FontSize', 14);
title('Avg. adjacent ciliary row spacing', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 


N = 19;
N2= 20;
data_ = table2array(readtable('WT-new.csv'));
midCount = [];

nCell = 1;
nCell2 = 1;
skipped = 0;
r = 2;
counts = zeros(300, N, 4);
centers = zeros(N, 1);
centers2 = zeros(N2, 1);

spacing = zeros(N, 3);
spacing2 = zeros(N2, 3);

spacing_ = zeros(338, 3);
nBBs = zeros(338, 1);
q = 1;
D = dir(main_path); 
for j=1:length(paths)
    full_path = fullfile(main_path, paths(j));
%         full_path = 'E:\BB Project\4- B1868_ _Poc1-mCh cell division';
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.

    for k = 3:length(D) % avoid using the first ones
        if D(k).isdir
            currD = D(k).name; % Get the current subdirectory name

            tempPath = fullfile(full_path, currD);

            
            data = csvread(fullfile(tempPath, 'intensityAnalysis.csv'));
            
            alignment = csvread(fullfile(tempPath, 'Alignment.csv'));
            spacing_data = csvread(fullfile(tempPath, 'Pairwise BB Distance Raw Data.csv'));
      
            spacing_data(isnan(spacing_data))=0;
            spacing_(q, 1) = mean(spacing_data(:, 1));
            spacing_(q, 2) = mean(spacing_data(:, 2));
            spacing_(q, 3) = mean(spacing_data(:, 3));
            nBBs(q) = size(data, 1);
            q = q + 1;
            
            if size(alignment, 1) == N
                nRow = size(alignment, 1);
                center_temp = zeros(N, 1);
                theta = data(:, 1);
                phi = data(:, 2);
    %             figure(8);

                for ii=1:nRow
                    idx_temp = [];
                    for jj=1:length(nonzeros(alignment(ii, :)))
                        if true
                            idx_temp(end+1) = alignment(ii, jj);
                        end
                    end
                    center_temp(ii) = std(data(idx_temp, 1));
                    if center_temp(ii) < 0
                        center_temp(ii) = 2*pi + center_temp(ii);
                    end
                    if center_temp(ii) > 2*pi
                        center_temp(ii) = center_temp(ii) - 2*pi;
                    end
    %                 scatter(sin(phi(idx_temp)).*cos(theta(idx_temp)), sin(phi(idx_temp)).*sin(theta(idx_temp)), 20 );
    %                 hold on
                end

                centers = centers + center_temp;
                nCell = nCell + 1;
                spacing = spacing + spacing_data;
            elseif size(alignment, 1) == N2
                nRow = size(alignment, 1);
                center_temp = zeros(N2, 1);
                theta = data(:, 1);
                phi = data(:, 2);
    %             figure(8);

                for ii=1:nRow
                    idx_temp = [];
                    for jj=1:length(nonzeros(alignment(ii, :)))

                        if true
                            idx_temp(end+1) = alignment(ii, jj);
                        end
                    end
                    center_temp(ii) = std(data(idx_temp, 1));
                    if center_temp(ii) < 0
                        center_temp(ii) = 2*pi + center_temp(ii);
                    end
                    if center_temp(ii) > 2*pi
                        center_temp(ii) = center_temp(ii) - 2*pi;
                    end
    %                 scatter(sin(phi(idx_temp)).*cos(theta(idx_temp)), sin(phi(idx_temp)).*sin(theta(idx_temp)), 20 );
    %                 hold on
                end

                centers2 = centers2 + center_temp;
                nCell2 = nCell2 + 1;
                spacing2 = spacing2 + spacing_data;
            end
        end
    end
end

centers = centers / nCell;
centers2 = centers2 / nCell2;
spacing = spacing / nCell;
spacing2 = spacing2 / nCell2;

spacing2(1, 1) = nan;
spacing2(20, 1) = nan;
spacing(1, 1) = nan;
spacing(19, 1) = nan;
figure();

% plot(centers, 'o-', 'LineWidth', 1.5);
% hold on
% plot(centers2, 'o-', 'LineWidth', 1.5);


plot((1:19)/19, spacing(:, [2, 3, 1]), 'o-', 'LineWidth', 1.5, 'HandleVisibility','off');
% yline(nanmean(spacing(:, 1)),'--','anterior avg','LineWidth',1.5, 'Color', '#0072BD');
% yline(nanmean(spacing(:, 2)),'--','posterior avg','LineWidth',1.5, 'Color', '#D95319');
% yline(nanmean(spacing(:, 3)),'--','medial avg','LineWidth',1.5, 'Color', '#EDB120');

% legend('anterior region','medial region', 'posterior region', 'FontSize', 14, 'Location', 'best');

xlabel('Row index', 'FontSize', 14);
ylabel('BB spacing (\mum)', 'FontSize', 14);
ax = gca; 
ax.ColorOrderIndex = 1;
hold on
plot((1:20)/20, spacing2(:, [2, 3, 1]), '*-', 'LineWidth', 1.5, 'HandleVisibility','off');
ax.ColorOrderIndex = 1;
hold on
plot((1:20)/20, spacing2(:, [2, 3, 1]), '-', 'LineWidth', 1.5);
% yline(nanmean(spacing2(:, 1)),'--','anterior avg','LineWidth',1.5, 'Color', '#0072BD');
% yline(nanmean(spacing2(:, 2)),'--','posterior avg','LineWidth',1.5, 'Color', '#D95319');
% yline(nanmean(spacing2(:, 3)),'--','medial avg','LineWidth',1.5, 'Color', '#EDB120');

ax.FontSize = 14; 
legend( 'posterior region', 'medial region', 'anterior region', 'FontSize', 14, 'Location', 'best');

xlabel('Normalized row index', 'FontSize', 14);

 % - Build title axes and title.
title('Avg. adjacent BB spacing', 'FontSize', 16);


figure();
params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
s = scatter(nBBs, spacing_(:, 2), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD','HandleVisibility','off');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;
hold on


y = lwppredict(nBBs, spacing_(:, 2), params, linspace(350, 720, 250)');
plot(linspace(350, 720, 250), y, 'LineWidth', 2, 'Color', '#0072BD');

hold on
s = scatter(nBBs, spacing_(:, 3), 'MarkerEdgeColor', '#D95319', 'MarkerFaceColor', '#D95319','HandleVisibility','off');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;
hold on


y = lwppredict(nBBs, spacing_(:, 3), params, linspace(350, 720, 250)');
plot(linspace(350, 720, 250), y, 'LineWidth', 2, 'Color', '#D95319');

hold on

y = lwppredict(nBBs, spacing_(:, 1), params, linspace(350, 720, 250)');
plot(linspace(350, 720, 250), y, 'LineWidth', 2, 'Color', '#EDB120');
hold on
s = scatter(nBBs, spacing_(:, 1), 'MarkerEdgeColor', '#EDB120', 'MarkerFaceColor', '#EDB120','HandleVisibility','off');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;



legend('posterior region', 'medial region', 'anterior region', 'FontSize', 14);
edges = 350:50:750;
set(gca, 'xtick', edges(1:end-1));

xlabel('Number of BBs', 'FontSize', 14);
ylabel('BB spacing (\mum)', 'FontSize', 14);
xlim([350 750]);
xticks(350:50:750);
ax = gca; 
ax.FontSize = 14
title('Avg. adjacent BB spacing through cell cycle', 'FontSize', 16);


