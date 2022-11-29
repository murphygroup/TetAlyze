clc
clear all


table = readtable('combined_data.csv', 'ReadVariableNames', true, 'ReadRowNames', true);
data = table2array(table);

t = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');

nexttile
h=histogram(data(:, 8), 25, 'Normalization','cdf')
y = h.Values';
x = h.BinEdges';
x(1:end-1) = x(1:end-1) + (x(2) - x(1))/2;
    
f=fit(y, x(1:end-1), 'smoothingspline', 'SmoothingParam',0.9995);    
d1 = differentiate(f, y);
h=figure();
yyaxis left
plot(f, y, x(1:end-1));
hold on
ylabel('Number of BBs', 'FontSize', 14);

yyaxis right
plot(y, d1);
xlabel('Cell cycle stage', 'FontSize', 14);
ylabel('Replication rate (derivative)', 'FontSize', 14);

ax = gca; 
ax.FontSize = 14; 