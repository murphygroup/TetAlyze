clc
close all

main_path = 'E:\BB data new\3. big1-1 Refeed experiment';   %  Path to the main results folder
hAll = [];
R2 = [];
Z = [];
w1All = [];
w2All = [];
labels = [];
numAll = [];
nRows = [];
nBBs = [];
X = [];
Y = [];
paths = ["N=1\B1868\0h"; "N=1\B1868\2h"; "N=1\B1868\4h"; "N=1\B1868\8h"; "N=1\B1868\24h";...
    "N=2\B1868\0h"; "N=2\B1868\2h"; "N=2\B1868\4h"; "N=2\B1868\8h"; "N=2\B1868\24h"; ...
     "N=3\B1868\0h"; "N=3\B1868\2h"; "N=3\B1868\4h"; "N=3\B1868\8h"; "N=3\B1868\24h";"4- B1868_ _Poc1-mCh cell division"; "temp"];
% paths = ["N=1\big1-1\0h"; "N=1\big1-1\2h"; "N=1\big1-1\4h"; "N=1\big1-1\8h"; "N=1\big1-1\24h";...
%     "N=2\big1-1\0h"; "N=2\big1-1\2h"; "N=2\big1-1\4h"; "N=2\big1-1\12h"; "N=2\big1-1\24h";...
%     "N=3\big1-1\0h"; "N=3\big1-1\2h"; "N=3\big1-1\4h"; "N=3\big1-1\12h"; "N=3\big1-1\24h";];
%     
% paths = ["N=1\B1868\8h", "N=2\B1868\8h", "N=3\B1868\8h"];
skipped = 0;
currLabel = 0;

n=0
for i=1:length(paths)
    full_path = fullfile(main_path, paths(i));
    D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.
    
    currLabel = currLabel + 1;
    if currLabel > 5
        currLabel = currLabel -5;
    end
    for k = 3:length(D) % avoid using the first ones
        currD = D(k).name; % Get the current subdirectory name
        n = n+1;
        tempPath = fullfile(full_path, currD);
        fList = dir(tempPath); % Get the file list in the subdirectory
        for j = 3:length(fList)
            currD = fList(j).name;
            if contains(currD, '.fig') && ~contains(currD, '._')
                [R2, Z, w1All, w2All, numAll, hAll, nRows, nBBs, X, Y] = readData(fullfile(tempPath, currD), R2, Z, w1All, w2All, numAll, hAll, nRows, nBBs, X, Y);
                labels = [labels currLabel];
            end
        end

    end
end
f=fit([Z' nBBs'], R2','poly32')
figure();
plot(f);
view([90 -90]);
xlabel({'$\bar{z}$'},'Interpreter','latex');
ylim([0 1.6]);
ylabel('# BBs');
zlabel('r^2');
% writematrix([w1All' w2All' numAll' hAll' labels' nRows'],'w-big1.csv'); 


f2=fit(Y', X','poly32')
figure();plot(f2);
zlabel({'$\bar{z}$'},'Interpreter','latex');
xlabel('empirical cummulative density');
ylabel('# BBs');

function [R2, Z, w1All, w2All, numAll, hAll, nRows, nBBs, X, Y] = readData(figname, R2, Z, w1All, w2All, numAll, hAll, nRows, nBBs, X, Y)
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
    data = vertcat(x,y)';
    [coeff,score,latent] = pca(data);
    nRows = [nRows length(h)];
    
    w1 = (max(score(1:end, 1))-min(score(1:end, 1)))/2;
    w2 = (max(score(1:end, 2))-min(score(1:end, 2)))/2;
    h = (max(z)-min(z))/2;
    a = score(1:end, 1)/w1;
    b = score(1:end, 2)/w2;
    r2 = (score(1:end, 1)/w1).^2 + (score(1:end, 2)/w2).^2;
    Z = [Z, (z-min(z))/(max(z)-min(z))];
    R2 = [R2 r2'];
    w1All = [w1All w1];
    w2All = [w2All w2];
    hAll = [hAll h];
    numAll = [numAll length(x')];
    n0 = length(x')/100;
    nBBs = [nBBs n0*ones(1, length(x'))];
    
        
    h=histogram((z-min(z))/(max(z)-min(z)), 150, 'Normalization', 'cdf');
    y = h.Values';
    x = h.BinEdges';
    x = x + x(2)/2;
    X = [X, x(1:end-1)'];
    Y = [Y, [y n0*ones(length(y), 1)]'];
end