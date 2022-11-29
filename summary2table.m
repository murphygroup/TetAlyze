clc
clear all
main_path = 'E:\BB Project\generative model images\2- Cell cycle generative model\N=3';   %  Path to the folder containing individual result folder
filename_ = 'N=3.csv';  % name of table file

data = [];
rowNames = {};
skipped = 0;

D = dir(main_path); 
for k = 3:length(D) % avoid using the first ones
    currD = D(k).name; % Get the current subdirectory name

    tempPath = fullfile(main_path, currD);
    if isfolder(tempPath)
        rowNames{end+1} = currD;
        fList = dir(tempPath); % Get the file list in the subdirectory
        for j = 3:length(fList)
            currD = fList(j).name;
            if contains(currD, '.fig') && ~contains(currD, '._')
                p = measurePerimeter(fullfile(tempPath, currD));
                break
            end
        end

        for j = 3:length(fList)
            currD = fList(j).name;
            if strcmp(currD, 'Summary.txt')
                T = readData(fullfile(tempPath, currD), p);
                data(end+1,:) = T;
            end
        end
    end
end
T_ = array2table(data, 'VariableNames', {'Anterior Area(um^2)', 'Anterior BB#', 'Medial Area(um^2)',...
    'Medial BB#', 'Posterior Area(um^2)', 'Posterior BB#', 'Cell surface area(um^2)', 'Assigned BB#', 'Ciliary row#', 'Avergae number of BBs per row',...
    'Cell height(um)', 'Cell major widths(um)', 'Cell minor widths(um)', 'Horizontal section perimeter(um)', 'Cell volume(um^3)', 'VolumeDeficit', 'Anterior neighbor BB pairwise distance(um)',...
    'Medial neighbor BB pairwise distance(um)', 'Posterior neighbor BB pairwise distance(um)', 'Anterior row pairwise distance(um)', 'Medial row pairwise distance(um)',...
    'Posterior row pairwise distance(um)', 'Mean of ciliary row lengths', 'STD of ciliary row lengths',  ...
     'Posterior immature BB#','Posterior medial immature BB#', 'Anterior medial immature BB#','Anterior immature BB#'},...
    'RowNames', rowNames);
writetable(T_, filename_, 'WriteVariableNames', true, 'WriteRowNames', true);


function T = readData(txtPath, p_)
q = 1;
T = zeros(1,28);

Str = fileread(txtPath);
Keys = ["Anterior Area (um^2):"; "# Anterior BBs:"; "Medial Area (um^2):"; "# Medial BBs:"; "Posterior Area (um^2):";...
    "# Posterior BBs:"; "# assigned BBs:"; "# BB rows:"; "Avergae number of BBs per row:"; "Cell height(um):"; "Cell widths(um):";...
    "Cell volume(um^3):"; "VolumeDeficit:"; "Average neighbor BB pairwise distance(um) (anterior, medial, posterior):"; "Average BB row pairwise distance(um) (anterior, medial, posterior):";...
    "Mean and standard deviation of ciliary row lengths:"; ];

    for p=1:length(Keys)
        Key = convertStringsToChars(Keys(p));
        Index = strfind(Str, Key);
        if p==11 
            Value = sscanf(Str(Index(1) + length(Key):end), '%g, %g', 2);
            for j=1:2
                T(q)=Value(j);
                q = q+1;
            end
            T(q)=p_;
            q = q+1;
        elseif  p==16
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
%         elseif p ==17
%             Value = sscanf(Str(Index(1) + length(Key):end), '%g, %g, %g, %g', 4);
%             for j=1:4
%                 T(q)=Value(j);
%                 q = q+1;
%             end
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
        if q == 7
            T(q) = T(1) + T(3) + T(5);
            q = q + 1;
            continue
        end
    end
end


function sum_ = measurePerimeter(figname)
    if (exist(figname,'file'))
        h1=openfig(figname, 'invisible');
    end
    ax = gca; 
    h = findobj(gca,'Type','line'); 

    x = [];
    y = [];
    z = [];

    for i=1:length(h)   
        x = [x h(i).XData];
        y = [y h(i).YData];
        z = [z h(i).ZData];
    end
    
    close(h1)
    data = vertcat(x,y)';
    [coeff,score,latent] = pca(data);
    k = convhull(score(:,1:2));
%     figure()
    sum_ = 0;
    for i=1:length(k)-1
        sum_ = sum_ + norm(score(k(i),1:2) - score(k(i+1),1:2));
%         plot(score(k(i:i+1),1),score(k(i:i+1),2));
%         hold on
    end
end
