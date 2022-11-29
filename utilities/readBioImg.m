%% [Mehul Bapat] (25 March 2020) Adding Relevant Documentation Below
    % data = bfopen("file_name") function returns an n-by-4 cell array
    %   data{s, 1} = m-by-2 cell array. (Represented by all_planes in the
    %   code)
    %   data{s, 2} = Original Metadata key/value pairs for s-th series
    %   data{s, 3} = Colour Look-Up tables for each plane in the s-th series
    %   data{s, 1}{t, 1} = Pixel Data for Plane t in series s.
    %   data{s, 1}{t, 2} = Label for plane t in the series s.
    % channel - The channel from which data is to be extracted
    % channel_count - The total number of channels present in the image.
    
%% [Mehul Bapat] (02/12/2020) Added vectorized method for extracting data
% Vectorized the whole I array/matrix without having to use nested for
% loops.
 %% [Ben] (12/06/17)
% Seems to be extracting pixel data from 3D image, and then reducing the
% dimensions of the representation by sampling every other plane and then
% reconstructing.

%% The Complete Function
function I=readBioImg(fileName, channel,channel_count)
%% First Reading the data from the file
data=bfopen(fileName); %This command opens the file. Full Path of the file recommended to be passed
% data = Tiff(fileName,'r');

%% Iterating for the pixel data over all the planes in all the series
all_planes = data{1, 1};% Cell array for all the planes in current series 1
total_np = size(all_planes, 1); % Total number of planes in the current plane dataset
np = total_np/channel_count; % Number of planes to be extracted from a particular channel

[m,n]=size(data{1,1}{1,1});
start_plane = channel;
end_plane = total_np;
req_planes = start_plane:channel_count:end_plane;
I=double(zeros(m,n,np));   % Preallocating the memory for better runtime.

% Vectorized Code for better Runtime
for i=1:length(req_planes)
    plane_temp = double(all_planes{req_planes(i), 1});%temporary variable for the current plane in the loop
%     plane_temp = plane_temp - min(plane_temp(:)) ;
%     plane_temp = plane_temp / max(plane_temp(:)) ;
    I(:,:,i) = plane_temp;
end

disp('Size of the output array = ')
disp(size(I));

end


















%% Old codes for new and old images
    %% Redundant Code (Vectorized code is written below)
%         [m,n]=size(plane_temp);% get size of the current plane
%         for i=1:m
%             for j=1:n
%                 I(i,j,floor((t+1)/2))=plane_temp(i,j);
%                 I(i,j,plane_count)=plane_temp(i,j);
%             end
%         end
    %% Old Code (Probably) for multiple-Channel Image Files
%channel_count = size(series1(2))
%series1_plane1=double(series1{1, 1}); % change type: uint16 to double
%[m,n]=size(series1_plane1) % get size of single plane
%series1_planeCount/2
% intensities in data_3d are the same as original dataset
%I=zeros(m,n,series1_planeCount/2); % only consider one channel
% I=zeros(m,n,series1_planeCount/2); % only consider one channel
% size(I)%for plane=channel:2:series1_planeCount
% for plane=channel:channel_count:series1_planeCount
%     plane
%     currPlane = double(series1{plane, 1});
%     series1{plane,2}
%     for i=1:m
%         for j=1:n
%             I(i,j,floor((plane+1)/2))=currPlane(i,j);
%         end
%     end
% end
