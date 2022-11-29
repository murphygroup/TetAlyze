clc
clear all
N = 19;
data_ = table2array(readtable('WT-new.csv'));
main_path = '2- Cell cycle generative model'; 
paths = ["N=1"; "N=2"; "N=3";];

h = figure
t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
ax1 = nexttile;
cmap = jet(20);
data19 = zeros(19, 250);
cmap = jet(20);
data20 = zeros(20, 250);

for i=1:19
    [newBBnum, rowBBnum, BBnum] = oneRow(main_path, paths, i, 19);
    params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
    data19(i, :) = lwppredict(BBnum', newBBnum'./rowBBnum', params, linspace(350, 600, 250)');
end

hp_all1 = zeros(19, 1);
for i=1:19
    if i==1
        nums = (data19(19, :) + data19(1, :) + data19(2, :))/3;
    elseif i==19
        nums = (data19(19, :) + data19(1, :) + data19(18, :))/3;
    else
        nums = (data19(i-1, :) + data19(i, :) + data19(i+1, :))/3;
    end
    plot(linspace(350, 600, 250)', nums, 'LineWidth', 2,  'Color', cmap(i, :));
    hp = find_half_point(nums);
    hp_all1(i) = hp;
    hold on
    plot([350 + hp, 350 + hp], [0, nums(hp)], '--', 'LineWidth', 1.5,  'Color', cmap(i, :));
    hold on
end
hold off
xlim([350, 600])
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Fraction of new BBs', 'FontSize', 14);
title('19 ciliary row cell population', 'FontSize', 14);
ax2 = nexttile;


for i=1:20
    [newBBnum, rowBBnum, BBnum] = oneRow(main_path, paths, i, 20);
    params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
    data20(i, :) = lwppredict(BBnum', newBBnum'./rowBBnum', params, linspace(350, 600, 250)');
end

hp_all2 = zeros(20, 1);


for i=1:20
    if i==1
        nums = (data20(20, :) + data20(1, :) + data20(2, :))/3;
    elseif i==20
        nums = (data20(20, :) + data20(1, :) + data20(19, :))/3;
    else
        nums = (data20(i-1, :) + data20(i, :) + data20(i+1, :))/3;
    end
    plot(linspace(350, 600, 250)', nums, 'LineWidth', 2,  'Color', cmap(i, :));
    hp = find_half_point(nums);
    hp_all2(i) = hp;
    hold on
    plot([350 + hp, 350 + hp], [0, nums(hp)], '--', 'LineWidth', 1.5,  'Color', cmap(i, :),'HandleVisibility','off');
    hold on
end
hold off

xlim([350, 600])
lgd = legend({'row 1', 'row 2', 'row 3', 'row 4', 'row 5', 'row 6', 'row 7', 'row 8', 'row 9', 'row 10',...
    'row 11', 'row 12', 'row 13', 'row 14', 'row 15', 'row 16', 'row 17', 'row 18', 'row 19' 'row 20'}, 'location','eastoutside','FontSize', 14);
lgd.NumColumns = 5;
xlabel('Number of BBs', 'FontSize', 14);
ylabel('Fraction of new BBs', 'FontSize', 14);
title('20 ciliary row cell population', 'FontSize', 14);
linkaxes([ax1, ax2], 'y');
ax1.FontSize = 14; 
ax2.FontSize = 14; 

lgd.Layout.Tile = 3;


figure()
plot((1:19)/19, 350+hp_all1, 'o-', 'LineWidth', 1.5);
hold on
plot((1:20)/20, 350+hp_all2, 'o-', 'LineWidth', 1.5);
xlabel('Normalized row index', 'FontSize', 14);
ylabel({'Number of BBs at ciliary rows'; 'with 50% of their BBs assembled'})
legend('19 ciliary rows','20 ciliary rows', 'FontSize', 14);
ax = gca; 
ax.FontSize = 14; 


table = readtable('combined_data.csv', 'ReadVariableNames', true, 'ReadRowNames', true);
data = table2array(table);


[y, x]=histcounts(data(:, 8), 25, 'Normalization','cdf');
y = 3 * y';
x = x';
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
xlabel('Relative time (fraction of cell cycle)', 'FontSize', 14);
ylabel('BB replication rate (BB/h)', 'FontSize', 14);

ax = gca; 
ax.FontSize = 14; 

function [newBBnum, rowBBnum, BBnum] = oneRow(main_path, paths, row_, N)
    r = 2;
    newBBnum = [];
    BBnum = [];
    rowBBnum = [];
    D = dir(main_path); 
    for j=1:length(paths)
        full_path = fullfile(main_path, paths(j));
        D = dir(full_path);

        for k = 3:length(D) 
            if D(k).isdir
                currD = D(k).name; 

                tempPath = fullfile(full_path, currD);
                data = csvread(fullfile(tempPath, 'intensityAnalysis.csv'));
                alignment = csvread(fullfile(tempPath, 'Alignment.csv'));
                if size(alignment, 1) ~= N
                    continue
                end
                nRow = size(alignment, 1);
                rowBBnum(end + 1) = sum(alignment(row_, :) > 0);
                
                idx_ =  find(data(:, 4)>r);
                countTemp = zeros(size(alignment, 1), 4);

                BBnum(end + 1) = size(data, 1);
                temp = 0;
                for ii=1:length(idx_)
                    [row ~]=find(alignment == idx_(ii));
                    if row == row_
                        temp =temp + 1;
                    end
                end
                newBBnum(end + 1) = temp;
            end
        end
    end
end


function hp = find_half_point(nums)
    total = sum(nums);
    running_sum = 0;
    for i=1:length(nums)
        running_sum = running_sum + nums(i);
        if running_sum > 0.5 * total
            hp = i;
            break
        end
    end
end