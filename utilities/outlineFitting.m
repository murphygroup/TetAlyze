clc
close all


outlineFuncs = {};
cdfFuncs = {};
nBBs = [];

hAll = [];
w1All = [];
w2All = [];
nRows = [];
stdAll = [];




main_path = 'E:\BB Project\generative model images\2- Cell cycle generative model';   %  Path to the main results folder

paths = ["N=1"; "N=2"; "N=3";];
skipped = 0;
currLabel = 0;

n=0


for i=1:length(paths)
    full_path = fullfile(main_path, paths(i));
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.

    for k = 3:length(D) % avoid using the first ones
        currD = D(k).name; % Get the current subdirectory name
        n = n+1;
        tempPath = fullfile(full_path, currD);
        fList = dir(tempPath); % Get the file list in the subdirectory
        for j = 3:length(fList)
            currD = fList(j).name;
            if contains(currD, '.fig') && ~contains(currD, '._')
                [outlineFuncs, cdfFuncs, nBBs, hAll, w1All, w2All, nRows, stdAll] = readData(fullfile(tempPath, currD), outlineFuncs, cdfFuncs, nBBs, hAll, w1All, w2All, nRows, stdAll);
            end
        end

    end
end
writematrix([nBBs' w1All' w2All' hAll' nRows' stdAll'],'WT-new.csv'); 
save('outlineFuncs.mat','outlineFuncs')
save('cdfFuncs.mat','cdfFuncs')


function p = epanechnikovKernel(x, mu, h)       
    p = (abs((x - mu)/h) <= 1).* (35/32 * (1-((x - mu)/h).^2).^3);
end
        
        
function [outlineFuncs, cdfFuncs, nBBs, hAll, w1All, w2All, nRows, stdAll] = readData(figname, outlineFuncs, cdfFuncs, nBBs, hAll, w1All, w2All, nRows, stdAll)
    if (exist(figname,'file'))
        h1=openfig(figname, 'invisible');
    end
    ax = gca; 
    h = findobj(gca,'Type','line'); 

    x = [];
    y = [];
    z = [];
    nRows = [nRows length(h)];
    BBnumPerRow = [];
    for i=1:length(h)   
        x = [x h(i).XData];
        y = [y h(i).YData];
        z = [z h(i).ZData];
        BBnumPerRow(end+1) = length(h(i).XData);
    end
    
    close(h1)
    std_ = std(BBnumPerRow);
    stdAll = [stdAll std_];
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
    
    w1All = [w1All w1];
    w2All = [w2All w2];
    hAll = [hAll h];
    
    
    h=histogram(z, 150, 'Normalization', 'cdf');
    y = h.Values';
    x = h.BinEdges';
    x = x + x(2)/2;
    
    f2=fit(y, x(1:end-1), 'poly4');
    cdfFuncs{end+1} = f2;

    % h=histogram((z-min(z))/(max(z)-min(z)), 150, 'Normalization', 'cdf');
    f=fit(z', r2,'smoothingspline', 'SmoothingParam', 0.995);
%     plot(f, z', r2);
    outlineFuncs{end+1} = f;
    nBBs(end+1) = length(score);
end
