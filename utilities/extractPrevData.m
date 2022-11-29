clc
close all


path = 'E:/BB Project/test'; % Path to a folder which contains several Alignment.fig files
% e.g.
% test
%   |- Alignment1.fig
%   |- Alignment2.fig
%    ...
resultPath = 'E:/BB Project/test/Summary.txt'  % Path to a text file, for writing results
D = dir(path);
fid = fopen(resultPath,'wt');


for k = 3:length(D) % avoid using the first ones
    currD = D(k).name; % Get the current subdirectory name
    if contains(currD, '.fig')
        tempPath = fullfile(path, currD);

        openfig(tempPath, 'invisible');
        ax = gca; 
        h = findobj(gca,'Type','line');

        cort_x = [];
        cort_y = [];
        cort_z = [];
        for i=1:length(h)
            l = h(i);
            for j=1:length(l.XData)
                cort_x(end+1) = l.XData(j);
                cort_y(end+1) = l.YData(j);
                cort_z(end+1) = l.ZData(j);
            end
        end
        antIdx = [];
        postIdx = [];
        midIdx = [];
        a = max(cort_z) ;
        b = min(cort_z);
        for i=1:length(cort_x)
            if cort_z(i) >= a - 0.25*(a-b)
                antIdx(end+1) = i;
            elseif cort_z(i) <= b + 0.25*(a-b)
                postIdx(end+1) = i;
            else
                midIdx(end+1) = i;
            end
        end
        cort_x = cort_x';
        cort_z = cort_z';
        cort_y = cort_y';
        P = alphaShape(cort_x,cort_y,cort_z);
        P.Alpha=20;
        figure(10)
        plot(P)
        antP = alphaShape(cort_x(antIdx),cort_y(antIdx),cort_z(antIdx));
        figure(11)
        plot(antP)
        antP.Alpha = 20;
        antArea = surfaceArea(antP);

        entireArea = surfaceArea(P);

        notAntIdx = setdiff(1:length(cort_x), antIdx);
        notAntP = alphaShape(cort_x(notAntIdx),cort_y(notAntIdx),cort_z(notAntIdx));
        notAntP.Alpha = 20;
        figure(12)
        plot(notAntP)
        notAntArea = surfaceArea(notAntP);
        tempArea1 = (notAntArea + antArea - entireArea)/2;


        postP = alphaShape(cort_x(postIdx),cort_y(postIdx),cort_z(postIdx));
        postP.Alpha = 20;
        figure(13)
        plot(postP)
        postArea = surfaceArea(postP);

        notPostIdx = setdiff(1:length(cort_x), postIdx);
        notPostP = alphaShape(cort_x(notPostIdx),cort_y(notPostIdx),cort_z(notPostIdx));
        notPostP.Alpha = 20;
        figure(14)
        plot(notPostP)

        notPostArea = surfaceArea(notPostP);
        tempArea2 = (notPostArea + postArea - entireArea)/2;

        antArea = antArea - tempArea1;
        postArea = postArea - tempArea2;

        midP = alphaShape(cort_x(midIdx),cort_y(midIdx),cort_z(midIdx));
        midP.Alpha = 20;
        midArea1 = surfaceArea(midP) - tempArea1 - tempArea2;

        fprintf(fid, '%s: \n', currD);
        fprintf(fid, 'Anterior Area (um^2): %f\n', antArea);
        fprintf(fid, '# Anterior BBs: %d\n', length(antIdx));

        fprintf(fid, 'Medial Area (um^2): %f\n', midArea1);
        fprintf(fid, '# Medial BBs: %d\n', length(midIdx));

        fprintf(fid, 'Posterior Area (um^2): %f\n', postArea);
        fprintf(fid, '# Posterior BBs: %d\n', length(postIdx));
        fprintf(fid, '\n');

        fprintf('%s: \n', currD);
        fprintf('Anterior Area (um^2): %f\n', antArea);
        fprintf('# Anterior BBs: %d\n', length(antIdx));

        fprintf('Medial Area (um^2): %f\n', midArea1);
        fprintf('# Medial BBs: %d\n', length(midIdx));

        fprintf('Posterior Area (um^2): %f\n', postArea);
        fprintf('# Posterior BBs: %d\n', length(postIdx));
        fprintf('\n');
    end
    
end