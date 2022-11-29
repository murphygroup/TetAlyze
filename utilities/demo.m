clc
close all
warning('off','all')

dir_path = 'E:\BB Project\generative model images\2. Cell cycle generative model\N=1';
fid = fopen(fullfile(dir_path, 'log.txt'),'w');


D = dir(dir_path);
for k = 3:length(D) % avoid using the first ones
    currD = D(k).name; % Get the current subdirectory name

    filename = fullfile(dir_path, currD);
    close all hidden
%     filename = 'E:\BB Project\generative model images\2. Cell cycle generative model\N=1\100522 B1868 Poc1-mCh asynchronous cells.sld - 13.tif' 


    %% parameters
    rejection_threshold = 2;  % 0.5 (strong rejection) - 3 (weak rejection) for WT or big1 without concaves, 3-7 for Big 1 with concaves
    minBBsInRow = 6; % Rows with fewer BBs than this value will be discarded.
    minRowLength = 5; % Threshold for BB row length.
    thresh_level = 2;

    %% June 2 2021, Paul Sun
    try
        I = readBioImg(filename, 1,1);
        imageID = '/';


        resultFolderPath = strrep(filename(1: end-4), '.', '-');

        mkdir(resultFolderPath)
        I2 = mat2gray(I);
        rng(1);
        [a, b, ~] = size(I2);
        a1 = (a - 650) / 2 + 1;
        a2 = a - a1;
        b1 = (b - 550) / 2 + 1;
        b2 = b - b1;
        vec = train_model(I2(a1:a2, b1:b2, :), thresh_level, rejection_threshold, false, minBBsInRow, minRowLength, resultFolderPath, imageID, I(a1:a2, b1:b2, :));

        f2 = join([resultFolderPath, '/data4PCA.csv']);
        writematrix(vec,f2);
    catch 
        fprintf(fid, filename(end-6:end-4));
        fprintf(fid, '\n');
    end
end
