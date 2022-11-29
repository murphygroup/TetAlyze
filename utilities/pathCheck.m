% path to the folder contains images and processed resutls
% e.g.
% /Users/lisamitchell/Documents/Pearson Lab Stuff/02052021 Poc1-mChy Cell Tag n=1/B1868 t=0
%   |--A.tif
%   |--A/...processed results
%   |--B.tif
%   |--B/...processed results
%   |--C.tif
%   |--C/...processed results
% path = '/Users/lisamitchell/Documents/Pearson Lab Stuff/02052021 Poc1-mChy Cell/Tag n=1/B1868 t=0'
path = 'E:\BB Project\B1868 t=0'
n = 0;
c = filesep;
D = dir(path);
for k = 3:length(D) % avoid using the first ones
    currD = D(k).name;
    if contains(currD, '.tif')
        tempPath = strrep(fullfile(path, currD(1:end-4)), '.', '-');
%         tempPath = strcat(tempPath, c);
        if ~isfolder(tempPath) && ~isfolder(fullfile(path, currD(1:end-4)))
            continue
        end
        n = n + 1;
    end
end

fprintf('%d processed result folders have been found.\n', n);