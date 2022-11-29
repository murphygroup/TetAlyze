clc
close all
clear all
rng(29);

% load data
data = table2array(readtable('combined_data.csv', 'ReadVariableNames', true, 'ReadRowNames', true));
load cdfFuncs.mat cdfFuncs
load outlineFuncs.mat outlineFuncs



for iii=1:15
    oneRun(iii, data, cdfFuncs, outlineFuncs, 19);
end
for iii=1:15
    oneRun(iii, data, cdfFuncs, outlineFuncs, 20);
end
function oneRun(id, data, cdfFuncs, outlineFuncs, nRow)
resultFolderPath = 'traj_models';
%% original

for i=1:nRow
    [newBBnum, rowBBnum, BBnum_] = countBBnum(i, nRow);
    if i == 1
        newBBnumAll = zeros(nRow, length(newBBnum));
        rowBBnumAll = zeros(nRow, length(rowBBnum));
    end
    newBBnumAll(i, :) = newBBnum;
    rowBBnumAll(i, :) = rowBBnum;
end
T1 = 5;
q = 0;
sig = false;
sig1= false;
nBB = max(370, normrnd(390, 5));
% nRow = 20;
nBBPR = (nBB + 14) / nRow;
nBBperRowFrac = [];
for i=1:nRow
    params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
    nBBperRowFrac(end + 1) = lwppredict(BBnum_', rowBBnumAll(i, :)', params, nBB);
end
nBBperRowFrac = [nBBperRowFrac(end), nBBperRowFrac];
nBBperRowFrac(end) = [];
nBBPR = (nBB + 14) * (nBBperRowFrac / sum(nBBperRowFrac));

BBNum = fix(normrnd(nBBPR', 1, [nRow 1]));

prob = normcdf(BBNum, nBBPR', 1);
BBNum(1) = BBNum(1) + 7;
BBNum(2) = BBNum(2) + 7;
N = sum(BBNum) - 14;
states = zeros(N, 1);
thetaAll = zeros(N, 1);
outAll = zeros(N, 1);

rAll = {};


theta = linspace(0, 2*pi-0.3, nRow);

xAll = [];
yAll = [];
zAll = [];

paramsW1 = lwpparams('EPA', 1, false, 0.2, 'robust', true);
w1 = 0.5 * lwppredict(data(:, 8), data(:, 12), paramsW1, N);

paramsW2 = lwpparams('EPA', 1, false, 0.2, 'robust', true);
w2 = 0.5 * lwppredict(data(:, 8), data(:, 13), paramsW2, N);

paramsH = lwpparams('EPA', 1, false, 0.2, 'robust', true);
h = lwppredict(data(:, 8), data(:, 11), paramsH, N);

alignment = zeros(nRow, max(BBNum));
thetaOld = [];
zOld = [];
nn = 0;
startPoints = max(0.03, norminv(1-prob, 0.03, 0.01));
endPoints = min(0.98, norminv(prob, 0.98, 0.01));


for l = 1:length(BBNum)
    n = BBNum(l);
    p = epanechnikovKernel(data(:, 8), N, 15);
    
    out = zeros(n, 1);

    for i=1:338
        f = cdfFuncs{i};
        f_return = f(linspace(startPoints(l), endPoints(l), n));
        f_return(isnan(f_return))=0;
        out = out + p(i) * f_return;
    end
    out = out/sum(p);

    out = out + normrnd(0, 0.025);

    out(out > 1) = 1;
    out(out < 0) = 0;

    r = zeros(n, 1);
    for i=1:338
        f = outlineFuncs{i};
        r =  r + p(i) * f(out);
    end
    r = sqrt(abs(r)/sum(p));

    r = r + normrnd(0, 0.025);

    thetaC = theta(l);
    x = [];
    y = [];
    z = [];

    tempTheta = [];
    tempZ = [];
    thetaLast = 0;

    for j = 1:n
        if l == 1
            if j < n - 6
                theta_ = normrnd(thetaC + (thetaLast-thetaC), 0.03);
                x(j) = real(w2 * r(j) * cos(theta_));
                y(j) = real(w1 * r(j) * sin(theta_));
                z(j) = out(j) * h;
                thetaLast = theta_;
                tempTheta(j) = theta_;
                tempZ(j) = z(j);
                nn = nn + 1;
                alignment(l, j) = nn;
                thetaAll(nn) = theta_;
                outAll(nn) = out(j);
            else
                tempTheta(j) = 0;
                tempZ(j) = out(j) * h;

            end
        elseif l == 2 && j > n - 7
            tempTheta(j) = 0;
            tempZ(j) = out(j) * h;
        else
            z(j) = out(j) * h;
            tempZ(j) = z(j);

            [~,I] = min(abs(z(j)-zOld));
            if j == 1
                thetaLast = thetaC;
            end
            theta_ = normrnd(thetaC + (thetaLast-thetaC) * 0.5 + 0.5 * thetaOld(I), 0.03);
            x(j) = real(w2 * r(j) * cos(theta_));
            y(j) = real(w1 * r(j) * sin(theta_));

            thetaLast = theta_;
            tempTheta(j) = theta_ - thetaC;
            nn = nn + 1;
            alignment(l, j) = nn;
            thetaAll(nn) = theta_;
            outAll(nn) = out(j);
        end
    end
    thetaOld = tempTheta(:,:);
    zOld = tempZ(:,:);
    xAll{l} = x;
    yAll{l} = y;
    zAll{l} = z;

end

% 
xAll_ = [];
yAll_ = [];
zAll_ = [];

for l=1:length(BBNum)

    xAll_ = [xAll_ xAll{l}];
    yAll_ = [yAll_ yAll{l}];
    zAll_ = [zAll_ zAll{l}];
end
filename = sprintf('%d_cycle_%d_%d.txt', nRow, id, q);
testAlignment(xAll_', yAll_', zAll_', alignment, filename, resultFolderPath);
q = q + 1;
% fig = figure();
for l = 1:length(BBNum)
    if l == 1 || l == 2
        n = BBNum(l)-7;
    else
        n = BBNum(l);
    end
%     plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12);
%     hold on
end
shp = alphaShape(0.95.*xAll_', 0.95.*yAll_', zAll_');
shp.Alpha = 15;
plot(shp, 'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.6,  'EdgeColor', [150 150 150]/255);
% 
axis equal
xlim([-10 10])
ylim([-18 18])
zlim([-1 55])
axis off

%% update
N0 = 720;
while N < N0
    n1 = 0;
    n2 = 0;
    n3 = 0;
    n4 = 0;
    
    for i=1:N
        if states(i) == 5
            states(i) = 0;
        elseif states(i) > 0
            states(i) = states(i) + 1;
        end
    end

    D = zeros(N, 1);
    all_ = [xAll_' yAll_' zAll_'];
    ii = 1;
    antIdx = [];
    medIdx1 = [];
    medIdx2 = [];
    postIdx = [];

    for l=1:size(alignment, 1)
        for i=1:length(find(alignment(l, :)))
            if i==1
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i+1), :), 'euclidean');
            elseif i == length(find(alignment(l, :)))
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i-1), :), 'euclidean');
            elseif sig1 && l == 1 && i == delIdx(1)
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i+1), :), 'euclidean');
            elseif sig1 && l == 1 && i == delIdx(1)-1
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i-1), :), 'euclidean');
            elseif sig && l == 2 && i == delIdx2(1)
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i+1), :), 'euclidean');
            elseif sig && l == 2 && i == delIdx2(1)-1
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i-1), :), 'euclidean');
            else
                D(alignment(l, i)) = pdist2(all_(alignment(l, i), :), all_(alignment(l, i+1), :), 'euclidean') * 0.5 + ...
                   0.5 * pdist2(all_(alignment(l, i), :), all_(alignment(l, i-1), :), 'euclidean');
            end

            if zAll_(alignment(l, i)) < 0.25 * h
                n4 = n4 + 1;
                if states(alignment(l, i)) == 0
                    postIdx(end+1) = alignment(l, i);
                end
            elseif zAll_(alignment(l, i)) < 0.5 * h
                n3 = n3 + 1;
                if states(alignment(l, i)) == 0
                    medIdx2(end+1) = alignment(l, i);
                end
            elseif zAll_(alignment(l, i)) < 0.75 * h
                n2 = n2 + 1;
                if states(alignment(l, i)) == 0
                    medIdx1(end+1) = alignment(l, i);
                end
            else
                n1 = n1 + 1;
                if states(ii) == 0
                    antIdx(end+1) = alignment(l, i);
                end
            end
            ii = ii + 1;
        end
    end

    rate1 = 0.0523;
    rate2 = 0.1639;
    rate3 = 0.1645;
    rate4 = 0.1361;

    n1 = round(rate1 * n1 / T1);
    n2 = round(rate2 * n2 / T1);
    n3 = round(rate3 * n3 / T1);
    n4 = round(rate4 * n4 / T1);


    for i=1:n4
        y = randsample(postIdx, 1, true, D(postIdx));
        [r c] = find(alignment == y);
        N = N + 1;

        if BBNum(r) == max(BBNum)
            alignment = [alignment zeros(nRow, 1)];
        end
        if c == BBNum(r) || (r == 1 && c == BBNum(r) - 7)|| (r == 2 && c == BBNum(r) - 7) ...
                || (sig && r == 1 && c == BBNum(r) - 7 - length(delIdx))|| (sig && r == 2 && c == BBNum(r) - 7 - length(delIdx2))
            thetaAll(end+1) = thetaAll(alignment(r, c)) + normrnd(0, 0.015);
        else
            thetaAll(end+1) = thetaAll(alignment(r, c)) * 0.5 + thetaAll(alignment(r, c + 1)) * 0.5 + normrnd(0, 0.015);
        end
        outAll(end+1) = min(1, outAll(alignment(r, c)) + 0.02);   
        BBNum(r) = BBNum(r) + 1;
        alignment(r, :) = [alignment(r, 1:c) N alignment(r,c+1:end-1)];

        postIdx(postIdx == y) = [];
        states(end+1) = 1;
        states(y) = 1;

        xAll_(end+1) = xAll_(alignment(r, c));
        yAll_(end+1) = xAll_(alignment(r, c));
        zAll_(end+1) = outAll(alignment(r, c)) * h;
        if sig && r == 1 && c < delIdx(1)
            delIdx = delIdx + 1;
        end
        if sig && r == 2 && c < delIdx(1)
            delIdx2 = delIdx2 + 1;
        end
    end

    for i=1:n3
        y = randsample(medIdx2, 1, true, D(medIdx2));
        [r c] = find(alignment == y);
        N = N + 1;

        if BBNum(r) == max(BBNum)
            alignment = [alignment zeros(nRow, 1)];
        end
        if c == BBNum(r) || (r == 1 && c == BBNum(r) - 7)|| (r == 2 && c == BBNum(r) - 7) ...
                || (sig && r == 1 && c == BBNum(r) - 7 - length(delIdx))|| (sig && r == 2 && c == BBNum(r) - 7 - length(delIdx2))
            thetaAll(end+1) = thetaAll(alignment(r, c)) + normrnd(0, 0.015);
        else
            thetaAll(end+1) = thetaAll(alignment(r, c)) * 0.5 + thetaAll(alignment(r, c + 1)) * 0.5 + normrnd(0, 0.015);
        end
        outAll(end+1) = min(1, outAll(alignment(r, c)) + 0.02);   
        BBNum(r) = BBNum(r) + 1;
        alignment(r, :) = [alignment(r, 1:c) N alignment(r,c+1:end-1)];

        medIdx2(medIdx2 == y) = [];
        states(end+1) = 1;
        states(y) = 1;

        xAll_(end+1) = xAll_(alignment(r, c));
        yAll_(end+1) = xAll_(alignment(r, c));
        zAll_(end+1) = outAll(alignment(r, c)) * h;
        if sig && r == 1 && c < delIdx(1)
            delIdx = delIdx + 1;
        end
        if sig && r == 2 && c < delIdx(1)
            delIdx2 = delIdx2 + 1;
        end
    end

    for i=1:n2
        y = randsample(medIdx1, 1, true, D(medIdx1));
        [r c] = find(alignment == y);
        N = N + 1;

        if BBNum(r) == max(BBNum)
            alignment = [alignment zeros(nRow, 1)];
        end
        if c == BBNum(r) || (r == 1 && c == BBNum(r) - 7)|| (r == 2 && c == BBNum(r) - 7) ...
                || (sig && r == 1 && c == BBNum(r) - 7 - length(delIdx))|| (sig && r == 2 && c == BBNum(r) - 7 - length(delIdx2))
            thetaAll(end+1) = thetaAll(alignment(r, c)) + normrnd(0, 0.015);
        else
            thetaAll(end+1) = thetaAll(alignment(r, c)) * 0.5 + thetaAll(alignment(r, c + 1)) * 0.5 + normrnd(0, 0.015);
        end
        outAll(end+1) = min(1, outAll(alignment(r, c)) + 0.02);   
        BBNum(r) = BBNum(r) + 1;
        alignment(r, :) = [alignment(r, 1:c) N alignment(r,c+1:end-1)];

        medIdx1(medIdx1 == y) = [];
        states(end+1) = 1;
        states(y) = 1;

        xAll_(end+1) = xAll_(alignment(r, c));
        yAll_(end+1) = xAll_(alignment(r, c));
        zAll_(end+1) = outAll(alignment(r, c)) * h;
        if sig && r == 1 && c < delIdx(1)
            delIdx = delIdx + 1;
        end
        if sig && r == 2 && c < delIdx2(1)
            delIdx2 = delIdx2 + 1;
        end
    end

    for i=1:n1
        y = randsample(antIdx, 1, true, D(antIdx));
        [r c] = find(alignment == y);
        N = N + 1;

        if BBNum(r) == max(BBNum)
            alignment = [alignment zeros(nRow, 1)];
        end
        if c == BBNum(r) || (r == 1 && c == BBNum(r) - 7)|| (r == 2 && c == BBNum(r) - 7) ...
                || (sig && r == 1 && c == BBNum(r) - 7 - length(delIdx))|| (sig && r == 2 && c == BBNum(r) - 7 - length(delIdx2))
            thetaAll(end+1) = thetaAll(alignment(r, c)) + normrnd(0, 0.015);
        else
            thetaAll(end+1) = thetaAll(alignment(r, c)) * 0.5 + thetaAll(alignment(r, c + 1)) * 0.5 + normrnd(0, 0.015);
        end
        outAll(end+1) = min(1, outAll(alignment(r, c)) + 0.02);   
        BBNum(r) = BBNum(r) + 1;
        alignment(r, :) = [alignment(r, 1:c) N alignment(r,c+1:end-1)];

        antIdx(antIdx == y) = [];
        states(end+1) = 1;
        states(y) = 1;

        xAll_(end+1) = xAll_(alignment(r, c));
        yAll_(end+1) = xAll_(alignment(r, c));
        zAll_(end+1) = outAll(alignment(r, c)) * h;
        if sig && r == 1 && c < delIdx(1)
            delIdx = delIdx + 1;
        end
        if sig && r == 2 && c < delIdx(1)
            delIdx2 = delIdx2 + 1;
        end
    end


    paramsW1 = lwpparams('EPA', 1, false, 0.20, 'robust', true);
    w1 = 0.5 * lwppredict(data(:, 8), data(:, 12), paramsW1, N);

    paramsW2 = lwpparams('EPA', 1, false, 0.20, 'robust', true);
    w2 = 0.5 * lwppredict(data(:, 8), data(:, 13), paramsW2, N);

    paramsH = lwpparams('EPA', 1, false, 0.20, 'robust', true);
    h = lwppredict(data(:, 8), data(:, 11), paramsH, N);

    p = epanechnikovKernel(data(:, 8), N, 20);
    
    for l = 1:length(BBNum)
        if l == 1
            if sig
                n = BBNum(l)-7-length(delIdx);
            else
                n = BBNum(l)-7;
            end
        elseif l == 2
            if sig
                n = BBNum(l)-7-length(delIdx2);
            else
                n = BBNum(l)-7;
            end
        else
            n = BBNum(l);
        end
        
        
        out = outAll(alignment(l, 1:n));
        
        if l == 1 && N > 650 && ~sig1
            delIdx = find(out>0.4 & out<0.5);
            out(delIdx) = [];
            alignment(l, delIdx(1):n-length(delIdx)) = alignment(l, delIdx(end)+1:n);
            alignment(l, n-length(delIdx)+1:n) = 0;
            n = n - length(delIdx);
            sig1 = true;
        end
        
        if l == 2 && N > 650 && ~sig
            delIdx2 = find(out>0.4 & out<0.5);
            out(delIdx2) = [];
            alignment(l, delIdx2(1):n-length(delIdx2)) = alignment(l, delIdx2(end)+1:n);
            alignment(l, n-length(delIdx2)+1:n) = 0;
            n = n - length(delIdx2);
            sig = true;

        end
        r = zeros(n, 1);
        
        
        zTarget = zeros(BBNum(l), 1);
        start_ = startPoints(l);
        end_ = endPoints(l);
        for i=1:338
            f = cdfFuncs{i};
            zTarget =  zTarget + p(i) * f(linspace(start_, end_, BBNum(l)));
        end
        zTarget = zTarget/sum(p) + normrnd(0, 0.025);
        if sig && l == 1 
            zTarget(delIdx) = [];
        end
        if sig && l == 2 
            zTarget(delIdx2) = [];
        end
        diff = (zTarget(1:n) - out)/T1;
        
        out = out + diff;

        for i=1:338
            f = outlineFuncs{i};
            r =  r + p(i) * f(out);
        end
        r = sqrt(r/sum(p));
        outAll(alignment(l, 1:n)) = out;
        

        r = r + normrnd(0, 0.025);
        for j=1:n
            xAll_(alignment(l, 1:n)) = real(w2 * r .* cos(thetaAll(alignment(l, 1:n))));
            yAll_(alignment(l, 1:n)) = real(w1 * r .* sin(thetaAll(alignment(l, 1:n))));
            zAll_(alignment(l, 1:n)) = out * h;
        end
    end
    filename = sprintf('%d_cycle_%d_%d.txt', nRow, id, q);
    testAlignment(xAll_', yAll_', zAll_', alignment, filename, resultFolderPath);
    q = q + 1;
    
end

% clf(fig,'reset')
    for l = 1:length(BBNum)
        if l == 1
            if sig
                n = BBNum(l)-7-length(delIdx);
    %             plot3(xAll_(alignment(l, 1:delIdx(1)-1)), yAll_(alignment(l, 1:delIdx(1)-1)), zAll_(alignment(l, 1:delIdx(1)-1)), '.-', 'MarkerSize',12);
    %             hold on
    %             plot3(xAll_(alignment(l, delIdx(1):n)), yAll_(alignment(l, delIdx(1):n)), zAll_(alignment(l, delIdx(1):n)), '.-', 'MarkerSize',12);
    %             hold on
            else
                n = BBNum(l)-7;
    %             plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12);
    %             hold on
            end
        elseif l == 2
            if sig
                n = BBNum(l)-7-length(delIdx2);
    %             plot3(xAll_(alignment(l, 1:delIdx2(1)-1)), yAll_(alignment(l, 1:delIdx2(1)-1)), zAll_(alignment(l, 1:delIdx2(1)-1)), '.-', 'MarkerSize',12);
    %             hold on
    %             plot3(xAll_(alignment(l, delIdx2(1):n)), yAll_(alignment(l, delIdx2(1):n)), zAll_(alignment(l, delIdx2(1):n)), '.-', 'MarkerSize',12);
    %             hold on
            else
                n = BBNum(l)-7;
    %             plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12);
    %             hold on
            end
        else
            n = BBNum(l);
    %         plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12);
    %         hold on
        end
    end

end
% shp = alphaShape(0.95.*xAll_', 0.95.*yAll_', zAll_');
% shp.Alpha = 15;
% plot(shp, 'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.6,  'EdgeColor', [150 150 150]/255);
% 
% axis equal
% xlim([-8 8])
% ylim([-16 16])
% zlim([-1 55])
% axis off






function p = epanechnikovKernel(x, mu, h)       
    p = (abs((x - mu)/h) <= 1).* (35/32 * (1-((x - mu)/h).^2).^3);
end