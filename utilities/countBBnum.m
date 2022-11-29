function  [newBBnum, rowBBnum, BBnum] = countBBnum(row_, N)
    main_path = '2- Cell cycle generative model'; 
    paths = ["N=1"; "N=2"; "N=3";];
    r = 2;
    newBBnum = [];
    BBnum = [];
    rowBBnum = [];
    D = dir(main_path); 
    for j=1:length(paths)
        full_path = fullfile(main_path, paths(j));
    %         full_path = 'E:\BB Project\4- B1868_ _Poc1-mCh cell division';
        D = dir(full_path); % A is a struct ... first elements are '.' and '..' used for navigation.

        for k = 3:length(D) % avoid using the first ones
            if D(k).isdir
                currD = D(k).name; % Get the current subdirectory name

                tempPath = fullfile(full_path, currD);


                data = csvread(fullfile(tempPath, 'intensityAnalysis.csv'));

                alignment = csvread(fullfile(tempPath, 'Alignment.csv'));
                if size(alignment, 1) ~= N
                    continue
                end
    %             dists = measurePairwiseDistance(fullfile(tempPath, 'Alignment.fig'));
                nRow = size(alignment, 1);
                rowBBnum(end + 1) = sum(alignment(row_, :) > 0);
                
                idx_ =  find(data(:, 4)>r);
                countTemp = zeros(size(alignment, 1), 4);

                BBnum(end + 1) = size(data, 1);
                temp = 0;
                for ii=1:length(idx_)
                    [row ~]=find(alignment == idx_(ii));
                    if row == row_
                        temp =temp + 1;
                    end
                end
                newBBnum(end + 1) = temp;



            end
        end
    end
end
