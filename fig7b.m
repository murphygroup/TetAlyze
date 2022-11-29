clc
clear all
main_path = '2- Cell cycle generative model';   %  Path to the main results folder

currLabel = 1;
labels = [];
data = [];
paths = ["N=1"; "N=2"; "N=3";];
names = ["measured"; "synthetic"];
skipped = 0;
for i=1:length(paths)
    full_path = fullfile(main_path, paths(i));
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.    
%     currLabel = currLabel + 1;
    for k = 3:length(D) % avoid using the first ones
        currD = D(k).name; % Get the current subdirectory name
        
        tempPath = fullfile(full_path, currD);
        fList = dir(tempPath); % Get the file list in the subdirectory
        for j = 3:length(fList)
            currD = fList(j).name;
            if strcmp(currD, 'Summary.txt')
%                 T = str2double(table2array(readtable(fullfile(tempPath, currD), 'ReadVariableNames', false)));
                T = readData(fullfile(tempPath, currD));
                if T(1) ~= 0
                    data(end+1,:) = T;
                    labels(end+1) = currLabel;
                else
                    skipped = skipped + 1;
                end
            end
        end

    end
end
% writematrix(data, 'data_all.csv');
% 
% full_path = 'E:\BB Project\generative modeling\models2';
% D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.    
% currLabel = currLabel + 1;
% for k = 3:length(D) % avoid using the first ones
%     currD = D(k).name; % Get the current subdirectory name
%     if contains(currD, '.txt')
%         T = readData(fullfile(full_path, currD));
%         if T ~= 0
%             data(end+1,:) = T;
%             labels(end+1) = currLabel;
%         else
%             skipped = skipped + 1;
%         end
%     end
% 
% end

lens = [15, 15, 15, 15, 15, 15, 15, 15, 15, 15];
full_path = 'traj_models';
D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.    
currLabel = currLabel + 1;
for k = 1:10 % avoid using the first ones
    for kk=1:lens(k)
        currD = sprintf('19_cycle_%d_%d.txt', k, kk-1);
        T = readData(fullfile(full_path, currD));

            data(end+1,:) = T;
            labels(end+1) = currLabel;

    end

end

lens2 = [17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17];
full_path = 'traj_models';
D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.    
currLabel = currLabel + 1;
for k = 1:10 % avoid using the first ones
    for kk=1:lens2(k)
        currD = sprintf('20_cycle_%d_%d.txt', k, kk-1);
        T = readData(fullfile(full_path, currD));

            data(end+1,:) = T;
            labels(end+1) = currLabel;

    end

end

% data2 = data(:, [1,3,5,7,13:end]);
for i = 1:22
    maxV = prctile(data(1:338, i), 95);
    minV = prctile(data(1:338, i), 5);
    


    data(1:338, i) = data(1:338, i) - mean(data(1:338, i));
    data(1:338, i) = data(1:338, i)/(maxV-minV);

    
    maxV = prctile(data(339:end, i), 95);
    minV = prctile(data(339:end, i), 5);


    data(339:end, i) = data(339:end, i) - mean(data(339:end, i));
    data(339:end, i) = data(339:end, i)/(maxV-minV);
end
[coeff,score,latent,tsquared,explained,mu] = pca(data(1:338, :));

score2 = data(339:end, :)*coeff;

latent2 = latent / sum(latent);
colors = {'green', 'magenta', 'cyan', 'red'};
markers = {'o'};
figure(10)

% color = colors{mod(i, 4)+1};
% marker = markers{fix(i/4) + 1};
scatter(score(1:338, 1), score(1:338, 2), 30,  colors{1}, 'o', 'HandleVisibility','off');
hold on

pointer = 1;
for k=1:10
    if k == 1
        plot(score2(pointer:pointer + lens(k)-1, 1), score2(pointer:pointer + lens(k)-1, 2), '-',  'Color', [1, 0, 1, 0.8]);
    else
        plot(score2(pointer:pointer + lens(k)-1, 1), score2(pointer:pointer + lens(k)-1, 2), '-',  'Color', [1, 0, 1, 0.8], 'HandleVisibility','off');
    end
    hold on
    pointer = pointer + lens(k);
end

for k=1:10
    if k == 1
        plot(score2(pointer:pointer + lens2(k)-1, 1), score2(pointer:pointer + lens2(k)-1, 2), '-', 'Color', [1, 0.6, 0., 0.8]);
    else
        plot(score2(pointer:pointer + lens2(k)-1, 1), score2(pointer:pointer + lens2(k)-1, 2), '-', 'Color', [1, 0.6, 0., 0.8], 'HandleVisibility','off');
    end
    hold on
    pointer = pointer + lens2(k);
end
% scatter(score3(:, 1), score3(:, 2), 30,  colors{2}, '*');
% hold on

% plot(score3(:, 1), score3(:, 2));
% scatter(score(labels==i+1, 1), score(labels==i+1, 2), 20,  'blue', 's');
% scatter(score(labels==i+2, 1), score(labels==i+2, 2), 20,  'red', 's');
% b = score(labels==i-1, 2)';
% title('PCA Embedding Plot')

xlabel(sprintf('PC1 (%.2f%%)', explained(1)), 'FontSize', 14);
ylabel(sprintf('PC2 (%.2f%%)', explained(2)), 'FontSize', 14);
box on
hold off

ax = gca; 
ax.FontSize = 14; 
legend('19 cilary rows','20 cilary rows', 'FontSize', 14);
warning('off','all')


function T = readData(txtPath)
q = 1;
T = zeros(1,19);

Str = fileread(txtPath);
Keys = ["Anterior Area (um^2):"; "# Anterior BBs:"; "Medial Area (um^2):"; "# Medial BBs:"; "Posterior Area (um^2):";...
    "# Posterior BBs:"; "# assigned BBs:"; "# BB rows:"; "Avergae number of BBs per row:"; "Cell height(um):"; "Cell widths(um):";...
    "Cell volume(um^3):"; "VolumeDeficit:"; "Average neighbor BB pairwise distance(um) (anterior, medial, posterior):"; "Average BB row pairwise distance(um) (anterior, medial, posterior):";...
    "Mean and standard deviation of ciliary row lengths:"];

    for p=1:length(Keys)
        Key = convertStringsToChars(Keys(p));
        Index = strfind(Str, Key);
        if p==11 || p==16
            Value = sscanf(Str(Index(1) + length(Key):end), '%g, %g', 2);
            for j=1:2
                T(q)=Value(j);
                q = q+1;
            end
        elseif p==14 || p==15
            Value = sscanf(Str(Index(1) + length(Key):end), '%g, %g, %g', 3);
            for j=1:3
                T(q)=Value(j);
                q = q+1;
            end
        else
            try
                Value = sscanf(Str(Index(1) + length(Key):end), '%f', 1);
            catch
                T = 0;
                return
            end
            T(q)=Value;
            q = q+1;
        end
    end
end




