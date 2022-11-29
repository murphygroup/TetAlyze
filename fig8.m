data = table2array(readtable('WT-new.csv'));

% params = lwpparams('EPA', 1, false);
% [hBest, critBest, results] = lwpfindh(data(:, 1), data(:, 2), params, 'CV');
params = lwpparams('EPA', 1, false, 0.2, 'robust', true);

figure
t = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(linspace(375, 720)', lwppredict(data(:, 1), data(:, 2), params, linspace(375, 720)'),'LineWidth',2);
hold on
plot(data(:, 1), data(:, 2), '.');
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Major axis widdth, w_1 (\mum)', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 
xticks(350:100:750);

nexttile
plot(linspace(375, 720)', lwppredict(data(:, 1), data(:, 3), params, linspace(375, 720)'),'LineWidth',2);
hold on
plot(data(:, 1), data(:, 3), '.');
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Minor axis width, w_2 (\mum)', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14;
xticks(350:100:750);

nexttile
plot(linspace(375, 720)', lwppredict(data(:, 1), data(:, 4), params, linspace(375, 720)'),'LineWidth',2);
hold on
plot(data(:, 1), data(:, 4), '.');
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Cell length, h (\mum)', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 
xticks(350:100:750);
nexttile


h1=openfig('E:\BB Project\generative model images\2- Cell cycle generative model\N=2\110522 B1868 Poc1-mCh asynchronous cells_1-sld - 4\Alignment.fig', 'invisible');

ax = gca; 
h = findobj(gca,'Type','line'); 

x = [];
y = [];
z = [];

BBnumPerRow = [];
for i=1:length(h)   
    x = [x h(i).XData];
    y = [y h(i).YData];
    z = [z h(i).ZData];
    BBnumPerRow(end+1) = length(h(i).XData);
end

close(h1)
std_ = std(BBnumPerRow);

data = vertcat(x,y)';
[coeff,score,latent] = pca(data);
% nRows = [nRows length(h)];

w1 = (max(score(1:end, 1))-min(score(1:end, 1)))/2;
w2 = (max(score(1:end, 2))-min(score(1:end, 2)))/2;
h = (max(z)-min(z))/2;
a = score(1:end, 1)/w1;
b = score(1:end, 2)/w2;
r2 = (score(1:end, 1)/w1).^2 + (score(1:end, 2)/w2).^2;
z = (z-min(z))/(max(z)-min(z));


h=histogram(z, 150, 'Normalization', 'cdf');
y = h.Values';
x = h.BinEdges';
x = x + x(2)/2;

% f2=fit(y, x(1:end-1), 'poly4');
% cdfFuncs{end+1} = f2;

% h=histogram((z-min(z))/(max(z)-min(z)), 150, 'Normalization', 'cdf');
f=fit(z', r2,'smoothingspline', 'SmoothingParam', 0.995);

plot(f,z', r2);
xlabel('Normalized cell length z''', 'FontSize', 14);
ylabel('r^2', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 