clear all
close all
table = readtable('N=1.csv', 'ReadVariableNames', true, 'ReadRowNames', true);
data1 = table2array(table);
table = readtable('N=2.csv', 'ReadVariableNames', true, 'ReadRowNames', true);
data2 = table2array(table);
table = readtable('N=3.csv', 'ReadVariableNames', true, 'ReadRowNames', true);
data3 = table2array(table);
table = readtable('combined_data.csv', 'ReadVariableNames', true, 'ReadRowNames', true);
data = table2array(table);

h=figure();
t = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');

nexttile
histogram(data(:, 8), 25, 'Normalization', 'probability');
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Fraction of cells', 'FontSize', 14);
xticks(350:50:750);
ax = gca; 
ax.FontSize = 14; 

nexttile
histogram(data(:, 9), 5, 'Normalization', 'probability');
xlabel('Number of ciliary rows', 'FontSize', 14);
ylabel('Fraction of cells', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 


nexttile
histogram(data(:, 1) + data(:, 3) + data(:, 5), 25, 'Normalization', 'probability');
xlabel('Cell surface area (\mum^2)', 'FontSize', 14);
ylabel('Fraction of cells', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 


nexttile
histogram(data(:, 15), 25, 'Normalization', 'probability');
xlabel('Cell volume (\mum^3)', 'FontSize', 14);
ylabel('Fraction of cells', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 

ax.XAxis.Exponent = 0
xtickformat('%.0f')
% 
% 
% h2=figure();
% t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
% nexttile
% scatter(data(:, 7), data(:, 1) + data(:, 3) + data(:, 5));
% xlabel('BB number', 'FontSize', 15);
% ylabel('Cell surface area (\mu^2)', 'FontSize', 15);
% r = corrcoef(data(:, 7), data(:, 1) + data(:, 3) + data(:, 5));
% text(350, 3500, sprintf('Pearson correlation = %0.3f', r(2, 1)), 'FontSize', 15);
% 
% nexttile
% scatter(data(:, 8), data(:, 13));
% xlabel('BB number', 'FontSize', 15);
% ylabel('Cell volume (\mu^3)', 'FontSize', 15);
% r = corrcoef(data(:, 8), data(:, 13));
% text(350, 15500, sprintf('Pearson correlation = %0.3f', r(2, 1)), 'FontSize', 15);
nexttile;

ax1 = nexttile;
% [~, edges] = histcounts(data(:, 8), 10);
% y = discretize(data(:, 8), edges);
% m = grpstats(data(:, 7), y, 'mean');
% % plot(edges(1:end-1), m, 'o-');
% std_ = grpstats(data(:, 7), y, 'std');

idx19 = find(data(:, 9) == 19);
idx20 = find(data(:, 9) == 20);
edges = 350:60:720;
% m = zeros(10, 1);
% std_ = zeros(10, 1);
% for i=1:10
%     idx = find(data(idx19, 8) >= edges(i) & data(idx19, 8) < edges(i+1));
%     m(i) = mean(data(idx, 7));
%     std_(i) = std(data(idx, 7));
% end
% 
% errorbar(edges(1:end-1)+20, m, std_);

s = scatter(data(:, 8), data(:, 7), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD','HandleVisibility','off');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;

hold on
params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
y = lwppredict(data(:, 8), data(:, 7), params, linspace(350, 660, 250)');
plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#0072BD');

hold on

% 
% s = scatter(data(idx20, 8), data(idx20, 7), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#D95319','HandleVisibility','off');
% s.AlphaData = 0.005;
% s.MarkerFaceAlpha = 0.2;
% s.MarkerEdgeAlpha = 0.2;
% hold on


% y = lwppredict(data(idx20, 8), data(idx20, 7), params, linspace(350, 660, 250)');
% plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#D95319');

set(gca, 'xtick', edges(1:end-1));
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Cell surface area (\mum^2)', 'FontSize', 14);
xlim([350 720]);
xticks(350:60:720);
ax = gca; 
ax.FontSize = 14; 


ax2 = nexttile;
% [~, edges] = histcounts(data(:, 8), 10);
% y = discretize(data(:, 8), edges);
% m = grpstats(data(:, 15), y, 'mean');
% % plot(edges(1:end-1), m, 'o-');
% std_ = grpstats(data(:, 15), y, 'std');
% errorbar(edges(1:end-1), m, std_);
s = scatter(data(:, 8), data(:, 15), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;

hold on
params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
y = lwppredict(data(:, 8), data(:, 15), params, linspace(350, 660, 250)');
plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#0072BD');

hold on

% 
% s = scatter(data(idx20, 8), data(idx20, 15), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#D95319');
% s.AlphaData = 0.005;
% s.MarkerFaceAlpha = 0.2;
% s.MarkerEdgeAlpha = 0.2;
% hold on
% 
% 
% y = lwppredict(data(idx20, 8), data(idx20, 15), params, linspace(350, 660, 250)');
% plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#D95319');
ax = gca;
ax.YRuler.Exponent = 0;
set(gca, 'xtick', edges(1:end-1));

xlabel('Number of BBs', 'FontSize', 14);
ylabel('Cell volume (\mum^3)', 'FontSize', 14);
xlim([350 720]);
xticks(350:60:720);
ax = gca; 
ax.FontSize = 14; 


nexttile;

idx19 = find(data(:, 9) == 19);
idx20 = find(data(:, 9) == 20);
edges = 350:30:720;

s = scatter(data(:, 8), data(:, 14), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD','HandleVisibility','off');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;

hold on
params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
y = lwppredict(data(:, 8), data(:, 14), params, linspace(350, 660, 250)');
plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#0072BD');

hold on


% s = scatter(data(idx20, 8), data(idx20, 14), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#D95319','HandleVisibility','off');
% s.AlphaData = 0.005;
% s.MarkerFaceAlpha = 0.2;
% s.MarkerEdgeAlpha = 0.2;
% hold on
% 
% 
% y = lwppredict(data(idx20, 8), data(idx20, 14), params, linspace(350, 660, 250)');
% plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#D95319');

set(gca, 'xtick', edges(1:end-1));
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Cell circumference (\mum)', 'FontSize', 14);
xlim([350 720]);
xticks(350:60:720);
ax = gca; 
ax.FontSize = 14; 


nexttile;
% [~, edges] = histcounts(data(:, 8), 10);
% y = discretize(data(:, 8), edges);
% m = grpstats(data(:, 15), y, 'mean');
% % plot(edges(1:end-1), m, 'o-');
% std_ = grpstats(data(:, 15), y, 'std');
% errorbar(edges(1:end-1), m, std_);
s = scatter(data(:, 8), data(:, 11), 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD','HandleVisibility','off');
s.AlphaData = 0.005;
s.MarkerFaceAlpha = 0.2;
s.MarkerEdgeAlpha = 0.2;

hold on
params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
y = lwppredict(data(:, 8), data(:, 11), params, linspace(350, 660, 250)');
plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#0072BD');

hold on


% s = scatter(data(idx20, 8), data(idx20, 11), 'MarkerEdgeColor', '#D95319', 'MarkerFaceColor', '#D95319','HandleVisibility','off');
% s.AlphaData = 0.005;
% s.MarkerFaceAlpha = 0.2;
% s.MarkerEdgeAlpha = 0.2;
% hold on
% 
% 
% y = lwppredict(data(idx20, 8), data(idx20, 11), params, linspace(350, 660, 250)');
% plot(linspace(350, 660, 250), y, 'LineWidth', 2, 'Color', '#D95319');
% legend('19 ciliary rows','20 ciliary rows', 'FontSize', 14, 'Location', 'best');
set(gca, 'xtick', edges(1:end-1));

xlabel('Number of BBs', 'FontSize', 14);
ylabel('Cell length (\mum)', 'FontSize', 14);
xlim([350 720]);
xticks(350:60:720);
ax = gca; 
ax.FontSize = 14; 



