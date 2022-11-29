clc
close all
clear all
rng(129);

% load data
data = table2array(readtable('combined_data.csv', 'ReadVariableNames', true, 'ReadRowNames', true));
load cdfFuncs.mat cdfFuncs
load outlineFuncs.mat outlineFuncs
filename = '19.gif';

nRow = 19;

nBBs = sort(data(:, 8));

%% original

% prepare OA pattern
del_idx = [53, 207, 70];
[im, map, alpha] = imread('OA.png');

im(:, :, 1) = im(:, :, 1).* uint8(alpha -50 > 0);
im(:, :, 2) = im(:, :, 2).* uint8(alpha -50 > 0) + 118 * uint8(alpha -50 == 0);
im(:, :, 3) = im(:, :, 3).* uint8(alpha -50 > 0) + 190 * uint8(alpha -50 == 0);
oa_pattern = im;
%----------------------------------------------------------------------------
for i=1:nRow
    [newBBnum, rowBBnum, BBnum_] = countBBnum(i, nRow);
    if i == 1
        newBBnumAll = zeros(nRow, length(newBBnum));
        rowBBnumAll = zeros(nRow, length(rowBBnum));
    end
    newBBnumAll(i, :) = newBBnum;
    rowBBnumAll(i, :) = rowBBnum;
end
T1 = 6;
q = 1;
sig = false;
sig1= false;
nBB =  362;
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
n0 = N;
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
startPoints = max(0.01, norminv(1-prob, 0.03, 0.01));
endPoints = min(0.98, norminv(prob, 0.98, 0.01));

OA = zeros(1, 3);
for l = 1:length(BBNum)
    n = BBNum(l);
    p = epanechnikovKernel(data(:, 8), N, 20);
    
    out = zeros(n, 1);
    
%     start_ = norminv(1-prob(l), 0.05, 0.02);
% %     start_ = normrnd(0.0, 0.05);
% %     end_ = normrnd(1, 0.02);
%     end_= norminv(prob(l), 0.99, 0.02);
%     startPoints(l) = start_;
%     endPoints(l) = end_;
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
        if ~ismember(i, del_idx)
            f = outlineFuncs{i};
            r =  r + p(i) * f(out);
        end
    end
    r = sqrt(abs(r)/sum(p));

    r = r + normrnd(0, 0.025);

    thetaC = theta(l);
    x = [];
    y = [];
    z = [];
    
    if l < 3
        OA = OA + [w2 * r(n - 3) * cos(theta(l)), w2 * r(3) * sin(theta(l)), out(n - 3) * h];
    end

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
% testAlignment(xAll_', yAll_', zAll_', alignment, q, id, nRow);
q = q + 1;
fig = figure();
for l = 1:length(BBNum)
    if l == 1 || l == 2
        n = BBNum(l)-7;
    else
        n = BBNum(l);
    end
    plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize', 12, 'Color', '#ff6700');
    hold on
end
OA = OA/2;
% plot3(OA(1), OA(2), OA(3), 'o', 'MarkerSize', 15, 'Color', '#1034A6', 'LineWidth', 3);
xImage = [OA(1) OA(1); OA(1)+1 OA(1)+1];   % The x data for the image corners
yImage = [OA(2)-3 OA(2)+3; OA(2)-3 OA(2)+3;];             % The y data for the image corners
zImage = [OA(3)+3 OA(3)+3; OA(3)-3 OA(3)-3;]; 
surf(xImage,yImage,zImage,...    % Plot the surface
     'CData',oa_pattern,...
     'FaceColor','texturemap','FaceAlpha','texturemap','LineStyle','none');
hold on 
shp = alphaShape(0.95.*xAll_', 0.95.*yAll_', zAll_');
shp.Alpha = 8;
plot(shp, 'FaceColor', '#0076BE', 'FaceAlpha', 0.7,  'EdgeColor', '#0076BE', 'EdgeAlpha', 0.7);
% plot3([3 3], [-17 -10], [1 1], 'o-', 'LineWidth', 2, 'Color', 'black');
hold on
% 
axis equal
% text(8, 16, 3, sprintf('BB#=%d',N));
% text(8, -10, 3, 'front side');
xlim([-10 10])
ylim([-18 18])
zlim([-1 55])
axis off
view([5 0 1])
% view([-5 0 1])
drawnow
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
% pause(5)

%% update
N0 = 720;
while  N < N0
    n1 = 0;
    n2 = 0;
    n3 = 0;
    n4 = 0;
    
    for i=1:N
        if states(i) == T1
            states(i) = 0;
        elseif states(i) > 0
            states(i) = states(i) + 1;
        end
    end

%     D = zeros(N, 1);
    all_ = [xAll_' yAll_' zAll_'];
    ii = 1;
    temp1 = {};
    temp2 = {};
    temp3 = {};
    temp4 = {};

    for l=1:size(alignment, 1)
        antIdx = [];
        medIdx1 = [];
        medIdx2 = [];
        postIdx = [];
        for i=1:length(find(alignment(l, :)))

            if zAll_(alignment(l, i)) < 0.25 * h
                if states(alignment(l, i)) == 0
                    postIdx(end+1) = alignment(l, i);
                end
            elseif zAll_(alignment(l, i)) < 0.5 * h
                if states(alignment(l, i)) == 0
                    medIdx2(end+1) = alignment(l, i);
                end
            elseif zAll_(alignment(l, i)) < 0.75 * h
                if states(alignment(l, i)) == 0
                    medIdx1(end+1) = alignment(l, i);
                end
            else
                if states(ii) == 0
                    antIdx(end+1) = alignment(l, i);
                end
            end
            ii = ii + 1;
        end
        temp4{end + 1} = postIdx;
        temp3{end + 1} = medIdx2;
        temp2{end + 1} = medIdx1;
        temp1{end + 1} = antIdx; 
    end

    rate1 = 0.0523;
    rate2 = 0.1639;
    rate3 = 0.1645;
    rate4 = 0.1361;

    n1 = round(rate1 * n1 / T1);
    n2 = round(rate2 * n2 / T1);
    n3 = round(rate3 * n3 / T1);
    n4 = round(rate4 * n4 / T1);
%     deltaN = nBBs(q) - nBBs(q-1);
    newBBperRowFrac = [];
    for i=1:nRow
        params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
        newBBperRowFrac(end + 1) = lwppredict(BBnum_', newBBnumAll(i, :)', params, N);
    end
    newBBperRowFrac = [newBBperRowFrac(end), newBBperRowFrac];
    newBBperRowFrac(end) = [];
    deltaN = sum(newBBperRowFrac)/T1; 
    newBBperRowFrac = newBBperRowFrac / sum(newBBperRowFrac);
    
    n1 = round(deltaN * rate1 / (rate1 + rate2 + rate3 + rate4));
    n2 = round(deltaN * rate2 / (rate1 + rate2 + rate3 + rate4));
    n3 = round(deltaN * rate3 / (rate1 + rate2 + rate3 + rate4));
    n4 = round(deltaN * rate4 / (rate1 + rate2 + rate3 + rate4));
    

    for i=1:n4
        r = randsample(1:nRow, 1, true, newBBperRowFrac');
        while isempty(temp4{r})
            r = randsample(1:nRow, 1, true, newBBperRowFrac');
        end
        y = randsample(temp4{r}, 1);
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
        r = randsample(1:nRow, 1, true, newBBperRowFrac');
        while isempty(temp3{r})
            r = randsample(1:nRow, 1, true, newBBperRowFrac');
        end
        y = randsample(temp3{r}, 1);
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
        r = randsample(1:nRow, 1, true, newBBperRowFrac');
        while isempty(temp2{r})
            r = randsample(1:nRow, 1, true, newBBperRowFrac');
        end
        y = randsample(temp2{r}, 1);
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
        r = randsample(3:nRow, 1, true, newBBperRowFrac(3:end)');
        while isempty(temp1{r})
            r = randsample(1:nRow, 1, true, newBBperRowFrac');
        end
        y = randsample(temp1{r}, 1);
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

    old_w1 = w1;
    old_w2 = w2;
    old_h = h;
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
        
        
    
        if l == 1 && N > 630 && ~sig1
            delIdx = find(out>0.4 & out<0.5);
            out(delIdx) = [];
            alignment(l, delIdx(1):n-length(delIdx)) = alignment(l, delIdx(end)+1:n);
            alignment(l, n-length(delIdx)+1:n) = 0;
            n = n - length(delIdx);
            sig1 = true;
        end
        
        if l == 2 && N > 630 && ~sig
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
            if ~ismember(i, del_idx)
                f = outlineFuncs{i};
                r =  r + p(i) * f(out);
            end
        end
        r = sqrt(r/(sum(p)-sum(p(del_idx))));
        outAll(alignment(l, 1:n)) = out;
        

        
        for j=1:n
            xAll_(alignment(l, 1:n)) = real(w2 * r .* cos(thetaAll(alignment(l, 1:n))));
            yAll_(alignment(l, 1:n)) = real(w1 * r .* sin(thetaAll(alignment(l, 1:n))));
            zAll_(alignment(l, 1:n)) = out * h;
        end
    end
%     testAlignment(xAll_', yAll_', zAll_', alignment, q, id, nRow);
    q = q + 1;
    


OA(1) = OA(1) * w1 / old_w1;
OA(2) = OA(2) * w2 / old_w2;
OA(3) = OA(3) * h / old_h;

clf(fig,'reset')
    for l = 1:length(BBNum)
        if l == 1
            if sig
                n = BBNum(l)-7-length(delIdx);
                plot3(xAll_(alignment(l, 1:delIdx(1)-1)), yAll_(alignment(l, 1:delIdx(1)-1)), zAll_(alignment(l, 1:delIdx(1)-1)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
                hold on
                plot3(xAll_(alignment(l, delIdx(1):n)), yAll_(alignment(l, delIdx(1):n)), zAll_(alignment(l, delIdx(1):n)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
                hold on
            else
                n = BBNum(l)-7;
                plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
                hold on
            end
        elseif l == 2
            if sig
                n = BBNum(l)-7-length(delIdx2);
                plot3(xAll_(alignment(l, 1:delIdx2(1)-1)), yAll_(alignment(l, 1:delIdx2(1)-1)), zAll_(alignment(l, 1:delIdx2(1)-1)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
                hold on
                plot3(xAll_(alignment(l, delIdx2(1):n)), yAll_(alignment(l, delIdx2(1):n)), zAll_(alignment(l, delIdx2(1):n)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
                hold on
            else
                n = BBNum(l)-7;
                plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
                hold on
            end
        else
            n = BBNum(l);
            plot3(xAll_(alignment(l, 1:n)), yAll_(alignment(l, 1:n)), zAll_(alignment(l, 1:n)), '.-', 'MarkerSize',12, 'Color', '#ff6700');
            hold on
        end 
    end
if exist('delIdx','var') == 1
    OA2 = [mean(xAll_([delIdx' delIdx2'])), mean(yAll_([delIdx' delIdx2'])), mean(zAll_([delIdx' delIdx2']))];
    xImage = [OA(1)+2 OA(1)+2; OA(1)+2.5 OA(1)+2.5];   % The x data for the image corners
    yImage = [OA(2)-3 OA(2)+3; OA(2)-3 OA(2)+3;];             % The y data for the image corners
    zImage = [0.45*h+3 0.45*h+3; 0.45*h-3 0.45*h-3;]; 
    surf(xImage,yImage,zImage,...    % Plot the surface
         'CData',oa_pattern,...
         'FaceColor','texturemap','FaceAlpha',0.8,'LineStyle','none');
    hold on 
end
xImage = [OA(1) OA(1); OA(1)+2 OA(1)+2];   % The x data for the image corners
yImage = [OA(2)-3 OA(2)+3; OA(2)-3 OA(2)+3;];             % The y data for the image corners
zImage = [OA(3)+3 OA(3)+3; OA(3)-3 OA(3)-3;]; 
surf(xImage,yImage,zImage,...    % Plot the surface
     'CData',oa_pattern,...
     'FaceColor','texturemap','FaceAlpha','texturemap','LineStyle','none');
hold on

shp = alphaShape(0.95.*xAll_', 0.95.*yAll_', zAll_');
shp.Alpha = 8;
plot(shp, 'FaceColor', '#0076BE', 'FaceAlpha', 0.7,  'EdgeColor', '#0076BE', 'EdgeAlpha', 0.7);
for l=1:N
    if states(l) == 1 || states(l) == 2
        color = '#ffd7b5';
    elseif states(l) == 3 || states(l) == 4
        color = '#ffb38a';
    elseif states(l) == 5
        color = '#ff9248';
    else
        continue
    end
    plot3(xAll_(l), yAll_(l), zAll_(l), '.-', 'MarkerSize',12, 'Color', color);
    hold on
end
% plot3([3 3], [-17 -10], [1 1], 'o-', 'LineWidth', 2, 'Color', 'black');
% hold on
% plot3([3], [-17 + sum(nBBs < N)/382 * 7], [1], '*', 'LineWidth', 2, 'Color', 'red');

axis equal
xlim([-8 8])
ylim([-18 18])
zlim([-1 55])
axis off
% text(8, 16, 3, sprintf('BB#=%d',N));
% text(8, -10, 3, 'front side');
view([5 0 1])
drawnow
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

imwrite(imind,cm,filename,'gif','WriteMode','append');

% pause(5);
end



function p = epanechnikovKernel(x, mu, h)       
    p = (abs((x - mu)/h) <= 1).* (35/32 * (1-((x - mu)/h).^2).^3);
end