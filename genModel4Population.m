clc
close all
rng(1);

% load data
data = table2array(readtable('combined_data.csv', 'ReadVariableNames', true, 'ReadRowNames', true));
load cdfFuncs.mat cdfFuncs
load outlineFuncs.mat outlineFuncs
resultFolderPath = '../single_models'
nRow = 19;
for i=1:nRow
    [newBBnum, rowBBnum, BBnum_19] = countBBnum(i, nRow);
    if i == 1
        rowBBnumAll_19 = zeros(nRow, length(rowBBnum));
    end
    rowBBnumAll_19(i, :) = rowBBnum;
end

nRow = 20;
for i=1:nRow
    [newBBnum, rowBBnum, BBnum_20] = countBBnum(i, nRow);
    if i == 1
        rowBBnumAll_20 = zeros(nRow, length(rowBBnum));
    end
    rowBBnumAll_20(i, :) = rowBBnum;
end

for ii=1:length(data)
    nBB = data(ii, 8) + 14;
    nRow = data(ii, 9);
    nBBPR = nBB / nRow;
   

    BBNum = fix(normrnd(nBBPR, 3, [nRow 1]));
    prob = normcdf(BBNum, nBBPR, 3);
    rAll = {};
    
    if nRow == 19
        nBBperRowFrac = [];
        for i=1:nRow
            params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
            nBBperRowFrac(end + 1) = lwppredict(BBnum_19', rowBBnumAll_19(i, :)', params, nBB);
        end
        nBBperRowFrac = [nBBperRowFrac(end), nBBperRowFrac];
        nBBperRowFrac(end) = [];
        nBBPR = (nBB + 14) * (nBBperRowFrac / sum(nBBperRowFrac));

        BBNum = fix(normrnd(nBBPR', 1, [nRow 1]));

        prob = normcdf(BBNum, nBBPR', 1);
        BBNum(1) = BBNum(1) + 7;
        BBNum(2) = BBNum(2) + 7;
    end
    
    if nRow == 20
        nBBperRowFrac = [];
        for i=1:nRow
            params = lwpparams('EPA', 1, false, 0.25, 'robust', true);
            nBBperRowFrac(end + 1) = lwppredict(BBnum_20', rowBBnumAll_20(i, :)', params, nBB);
        end
        nBBperRowFrac = [nBBperRowFrac(end), nBBperRowFrac];
        nBBperRowFrac(end) = [];
        nBBPR = (nBB + 14) * (nBBperRowFrac / sum(nBBperRowFrac));

        BBNum = fix(normrnd(nBBPR', 1, [nRow 1]));

        prob = normcdf(BBNum, nBBPR', 1);
        BBNum(1) = BBNum(1) + 7;
        BBNum(2) = BBNum(2) + 7;
    end

    N = sum(BBNum) - 14;
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
    temp = max(endPoints - startPoints);
    h = h / temp;
    for l = 1:length(BBNum)
        n = BBNum(l);
        p = epanechnikovKernel(data(:, 8), nBB-14, 20);
        out = zeros(n, 1);
%         start_ = endPoints;
%         end_ = min(1, normrnd(0.95, 0.03));
        for i=1:338
            f = cdfFuncs{i};
            out =  out + p(i) * f(linspace(startPoints(l), endPoints(l), n));
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
        r = sqrt(r/sum(p));

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
                if j < n - 7
                    theta_ = normrnd(thetaC + (thetaLast-thetaC), 0.03);
                    x(j) = normrnd(real(w2 * r(j) * cos(theta_)), 0.2);
                    y(j) = normrnd(real(w1 * r(j) * sin(theta_)), 0.2);
                    z(j) = out(j) * h;
                    thetaLast = theta_;
                    tempTheta(j) = theta_;
                    tempZ(j) = z(j);
                    nn = nn + 1;
                    alignment(l, j) = nn;
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
            end
        end
        thetaOld = tempTheta(:,:);
        zOld = tempZ(:,:);
        xAll{l} = x;
        yAll{l} = y;
        zAll{l} = z;

    end

% 
%     figure();
    xAll_ = [];
    yAll_ = [];
    zAll_ = [];
    for l=1:length(BBNum)
%         plot3(xAll{l}, yAll{l}, zAll{l}, 'o-');
%         hold on
        xAll_ = [xAll_ xAll{l}];
        yAll_ = [yAll_ yAll{l}];
        zAll_ = [zAll_ zAll{l}];
    end
    
%     shp = alphaShape(0.95.*xAll_', 0.95.*yAll_', zAll_');
%     shp.Alpha = 15;
%     plot(shp, 'FaceColor', [150 150 150]/255, 'FaceAlpha', 0.6,  'EdgeColor', [150 150 150]/255);
%     axis equal
    filename = sprintf('model_%d.txt', ii);
    alignment = testAlignment(xAll_', yAll_', zAll_', alignment, filename, resultFolderPath);

end

function p = epanechnikovKernel(x, mu, h)       
    p = (abs((x - mu)/h) <= 1).* (35/32 * (1-((x - mu)/h).^2).^3);
end