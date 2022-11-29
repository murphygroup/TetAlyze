clc
close all
warning('off','all')


filename = 'N=3/200522 B1868 Poc1-mCh asynchronous cells_1.sld - 12.tif'


rejection_threshold = 2.;  % 0.5 (strong rejection) - 3 (weak rejection) 
minBBsInRow = 5; % Rows with fewer BBs than this value will be discarded.
minRowLength = 8; % Threshold for BB row length.
thresh_level = 2;


I = readBioImg(filename, 1,1);
imageID = '/';


resultFolderPath = strrep(filename(1: end-4), '.', '-');

mkdir(resultFolderPath)
I2 = mat2gray(I);
rng(1);
[a, b, ~] = size(I2);
a1 = (a - 500) / 2 + 1;
a2 = a - a1;
b1 = (b - 500) / 2 + 1;
b2 = b - b1;
vec = train_model(I2(a1:a2, b1:b2, :), thresh_level, rejection_threshold, true, minBBsInRow, minRowLength, resultFolderPath, imageID, I(a1:a2, b1:b2, :));

f2 = join([resultFolderPath, '/data4PCA.csv']);
writematrix(vec,f2);
% close all
