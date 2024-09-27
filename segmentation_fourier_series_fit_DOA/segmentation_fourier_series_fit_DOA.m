clc
close all
clear all

%% Fourier series settings

total_number_of_series_cell = 15;

total_number_of_series_nuc = 5;

%% import images

% image resolution
x_res = 1/2.0292; %um/pixel
y_res = x_res; %um/pixel

% import cytoplasm segmentation 
I = imread('cytoplasm_segmentation.png');
I = I(2:end-1,2:end-1,:); % removes black border

% import nuclear segmentation 
DAPI = imread('nuclear_segmentation.png');
DAPI = DAPI(2:end-1,2:end-1,:); % removes black border

% import aSMA image
aSMA = imread('Alexa 488.tiff');    
aSMA = aSMA(2:end-1,2:end-1,:); % removes black border

% import manual activation gradine
activation = imread('activation_segmentation.png');
activation = activation(2:end-1,2:end-1,:); % removes black border
activation_r = activation(:,:,1);
activation_g = activation(:,:,2);
activation_b = activation(:,:,3);

% compute orientation and coherence
[orientation,coherence] = aSMA_orientation(aSMA,4);

% rescaling
orientation_t = rescale(orientation);
coherence_t = coherence;
aSMA_t = rescale(double(aSMA));

% generate orientation image
hsv_img = hsv2rgb(cat(3,orientation_t,coherence_t,aSMA_t));

% shows orientation images. Uncomment to see. 
% figure 
% imshow(hsv_img)

% creates plasma colormap
[cmap] = parula(256);

% Map the grayscale image to an indexed image based on the colormap size
indexedImage = gray2ind(coherence_t, size(cmap, 1));

% Convert the indexed image to RGB using the plasma colormap
rgbImage = ind2rgb(indexedImage, cmap);

% shows orientation images. Uncomment to see. 
% figure
% imshow(rgbImage)

% saves orientation and coherence images.
imwrite(rgbImage,"coherence.png");
imwrite(hsv_img,"orientation.png");

%% pre-allocating morphological variables and activation variable and DOA

cell_area = [];
cell_roundness = [];
cell_width = [];
cell_length = [];
cell_aspect_ratio = [];
cell_solidity = [];

DAPI_area = [];
DAPI_roundness = [];
DAPI_width = [];
DAPI_length = [];
DAPI_aspect_ratio = [];
DAPI_solidity = [];

activated = [];
DOA = []; % degree of anisotropy

%% segment individual cytoplasms

% concatenates red, green, and blue channels 
I_1 = I(:,:,1);
I_2 = I(:,:,2);
I_3 = I(:,:,3);

% concatenates red, green, and blue values into a list
I2 = [I_1(:),I_2(:),I_3(:)];

% look for unique pixel colors
C = unique(I2,'rows');

% remove black pixels which are sorted to the first row
C(1,:) = [];
 
% cell counter
cell_num = 1;

% pre-allocate cell to contain the pixel indices of every cell identified
% in the segmentation image
cell_indx = {};

% iterates for every unique cell color found in the segmentation images
for i = 1:length(C)

    % finds pixels that are the same color as C(i,:)
    indx1 = find(I_1 == C(i,1));
    indx2 = find(I_2 == C(i,2));
    indx3 = find(I_3 == C(i,3));
    indx_all = intersect(indx1,indx2);
    indx_all = intersect(indx_all,indx3);

    % extracts and binarizes cell body that is colored as C(i,:)
    bw = zeros(size(I_1));
    bw(indx_all) = 1;
    bw = bw > 0;
    bw = imclearborder(bw); % gets rid of cells touching border
    
    % if cells touch border, move on to next iteration
    if sum(sum(sum(bw))) < 1
        continue
    end
    
    % generate pixel indices for cells colored as C(i,:)
    stats = regionprops('table',bw,'PixelIdxList');

    % saves the pixel indices of the identified cell
    for j = 1:size(stats,1)
        cell_indx{cell_num} = stats.PixelIdxList{j,1};
        cell_num = cell_num + 1;
    end

end
    
%% segment individual nuclei
    
% concatenates red, green, and blue channels 
DAPI_1 = DAPI(:,:,1);
DAPI_2 = DAPI(:,:,2);
DAPI_3 = DAPI(:,:,3);

% concatenates red, green, and blue values into a list
DAPI2 = [DAPI_1(:),DAPI_2(:),DAPI_3(:)];

% look for unique pixel colors
C2 = unique(DAPI2,'rows');

% remove black pixels which are sorted to the first row
C2(1,:) = [];

% cell counter
cell_num = 1;

% pre-allocate cell to contain the pixel indices of every nuclei identified
% in the segmentation image
DAPI_indx = {};

% iterates for every unique nuclei color found in the segmentation images
for i = 1:length(C2)

    % finds pixels that are the same color as C2(i,:)
    indx1 = find(DAPI_1 == C2(i,1));
    indx2 = find(DAPI_2 == C2(i,2));
    indx3 = find(DAPI_3 == C2(i,3));
    indx_all = intersect(indx1,indx2);
    indx_all = intersect(indx_all,indx3);

    % extracts and binarizes nuclei that is colored as C2(i,:)
    bw = zeros(size(DAPI_1));
    bw(indx_all) = 1;
    bw = bw > 0;
    bw = imclearborder(bw); % gets rid of nuclei touching border
    
    % if nuclei touches border, move on to next iteration
    if sum(sum(sum(bw))) < 1
        continue
    end
   
    % generate pixel indices for nuclei colored as C(i,:)
    stats = regionprops('table',bw,'PixelIdxList');

    % saves the pixel indices of the identified nuclei
    for j = 1:size(stats,1)
        DAPI_indx{cell_num} = stats.PixelIdxList{j,1};
        cell_num = cell_num + 1;
    end

end
    
%% match nuclei and cytoplasm, compute morphological measures, shape registration, Fourier series fit

% cell counter
cell_num = 1;
    
% iterates for every unique cell identified 
for i = 1:length(cell_indx)

    % pre-allocate variables to store cytoplasm and dapi images
    CYTO_bw = zeros(size(I_1));
    DAPI_bw = zeros(size(DAPI_1));

    % extracts and binarizes current cell and computes morphology metrics
    cytoplasm = cell_indx{1,i};
    CYTO_bw(cytoplasm) = 1;
    CYTO_bw = CYTO_bw > 0;
    CYTO_stats = regionprops('table',CYTO_bw,'PixelIdxList','Area','Circularity','MajorAxisLength','MinorAxisLength','Solidity');

    % iterates through the nuclei indices to match the nuclei that overlaps
    % the most to the current cell 
    number = [];
    for j = 1:length(DAPI_indx)
        DAPI = DAPI_indx{1,j};
        [C2,ia,ib] = intersect(cytoplasm,DAPI);
        number(j) = length(ia);
    end
    
    % if no nuclei overlaps with current cell, move on to the next
    % iteration
    if isempty(nonzeros(number))
        continue
    end

    % finds the nuclei the overlaps the most with the current cell
    DAPI_num = find(number == max(number));

    % extracts and binarizes matched nuclei and computes morphology metrics
    DAPI_bw(DAPI_indx{1,DAPI_num}) = 1;
    DAPI_bw = DAPI_bw > 0;
    DAPI_stats = regionprops('table',DAPI_bw,'PixelIdxList','Area','Circularity','MajorAxisLength','MinorAxisLength','Solidity');

    % compile morphology metrices for current cell body
    cell_area = [cell_area;CYTO_stats.Area*x_res*y_res];
    cell_roundness = [cell_roundness;CYTO_stats.Circularity];
    cell_width = [cell_width;CYTO_stats.MinorAxisLength*x_res];
    cell_length = [cell_length;CYTO_stats.MajorAxisLength*x_res];
    cell_aspect_ratio = [cell_aspect_ratio;CYTO_stats.MajorAxisLength/CYTO_stats.MinorAxisLength];
    cell_solidity = [cell_solidity;CYTO_stats.Solidity];

    % compile morphology metrices for matched nuclei
    DAPI_area = [DAPI_area;DAPI_stats.Area*x_res*y_res];
    DAPI_roundness = [DAPI_roundness;DAPI_stats.Circularity];
    DAPI_width = [DAPI_width;DAPI_stats.MinorAxisLength*x_res];
    DAPI_length = [DAPI_length;DAPI_stats.MajorAxisLength*x_res];
    DAPI_aspect_ratio = [DAPI_aspect_ratio;DAPI_stats.MajorAxisLength/DAPI_stats.MinorAxisLength];
    DAPI_solidity = [DAPI_solidity;DAPI_stats.Solidity];

    % checks cell activation, assigns 0 for unactivated and 1 for activated
    activation_check = mean(activation_r(CYTO_bw))/mean([activation_g(CYTO_bw);activation_b(CYTO_bw)]);
    if activation_check > 1 % checks if cell is red (unactivated)
        activated = [activated;0];
    else
        activated = [activated;1];
    end

    % computes DOA for every identified cell
    DOA = [DOA;median(coherence_t(CYTO_bw))];

    %% Shape registration and Fourier fit for cell body
   
    % extracts boundary point of current cell
    B = bwboundaries(CYTO_bw);
    x = B{1,1}(:,1);
    y = B{1,1}(:,2);

    % removes redundant boundary points
    points = [x,y];
    points = unique(points,'rows','stable');

    % transforms boundary points from pixel indices to cartesian
    % coordinates
    points = [points(:,2),(size(CYTO_bw,1)-points(:,1))];

    % centers points
    points = points - mean(points);

    % aligns cell boundary so that it's longest axis is on the x-axis 
    R = pca(points);
    points = points*R;

    % converts points from pixels to um
    points_um = points*x_res;

    % rotates shape so that largest area is on the right of the y-axis
    positive = find(points_um(:,1)>0);
    negative = find(points_um(:,1)<0);
    left = max(points_um(negative,2)) - min(points_um(negative,2));
    right = max(points_um(positive,2)) - min(points_um(positive,2));

    if right > left
        theta = 0;
    else
        theta = 180;
    end

    R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];

    points_um = (R*points_um')';
    
    % flips shape horizontally if most of the boundary points are below x-axis
    top = find(points_um(:,2)>0);
    bot = find(points_um(:,2)<0);

    if length(bot) > length(top)
        points_um(:,2) = -points_um(:,2);
    end

    %controls winding to be counter-clockwise
    tf = ispolycw(points_um(:,1),points_um(:,2));
    if tf 
        points_um = flip(points_um,1);
    end

    % re-orders points such that point furthest in the x-axis is first on
    % the list of points
    indx_ = find(points_um(:,1) == max(points_um(:,1)));
    points_um = [points_um(indx_:end,:);points_um(1:indx_-1,:)];

    % finds arc-length of distance around cell shape
    dt = diff(points_um);
    dis = zeros(length(dt),1);
    for k = 1:length(dt)
        if k == 1
            dis(k) = norm(dt(k,:));
        else
            dis(k) = dis(k-1) + norm(dt(k,:));
        end
    end
    
    % prepares arc-length data for Fourier series fit
    dis = [0;dis]; % first point is zero arc-length
    dis = normalize(dis,'range')*2*pi; % normalize from arclength to 2pi
    
    % concatenates x and y points separately 
    points_um_x_CELL = points_um(:,1); % concatenates x
    points_um_y_CELL = points_um(:,2); % concatenates y
    
    % computes Fourier series coefficients
    [a] = fourier_series_fit(dis,points_um_x_CELL,total_number_of_series_cell);
    [b] = fourier_series_fit(dis,points_um_y_CELL,total_number_of_series_cell);

    % saves Fourier series coefficients
    save(strcat('a_coeff_CELL_',num2str(cell_num),'.mat'),'a')
    save(strcat('b_coeff_CELL_',num2str(cell_num),'.mat'),'b')

    % Fourier series evaluated at these theta values
    thetas = linspace(0,2*pi,10000);
    
    % Evaluation of Fourier series fit
    [x_predicted_CELL] = fourier_series_evaluate(a,dis);
    [y_predicted_CELL] = fourier_series_evaluate(b,dis);

    %% Shape registration and Fourier fit for nuclei
   
    % extracts boundary point of current nuclei
    B = bwboundaries(DAPI_bw);
    x = B{1,1}(:,1);
    y = B{1,1}(:,2);

    % removes redundant boundary points
    points = [x,y];
    points = unique(points,'rows','stable');

    % transforms boundary points from pixel indices to cartesian
    % coordinates
    points = [points(:,2),(size(DAPI_bw,1)-points(:,1))];

    % centers points
    points = points - mean(points);

    % aligns nuclei boundary so that it's longest axis is on the x-axis 
    R = pca(points);
    points = points*R;

    % converts points from pixels to um
    points_um = points*x_res;

    % rotates shape so that largest area is on the right of the y-axis
    positive = find(points_um(:,1)>0);
    negative = find(points_um(:,1)<0);
    left = max(points_um(negative,2)) - min(points_um(negative,2));
    right = max(points_um(positive,2)) - min(points_um(positive,2));

    if right > left
        theta = 0;
    else
        theta = 180;
    end

    R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];

    points_um = (R*points_um')';
    
    % flips shape horizontally if most of the boundary points are below x-axis
    top = find(points_um(:,2)>0);
    bot = find(points_um(:,2)<0);

    if length(bot) > length(top)
        points_um(:,2) = -points_um(:,2);
    end

    %controls winding to be counter-clockwise
    tf = ispolycw(points_um(:,1),points_um(:,2));
    if tf 
        points_um = flip(points_um,1);
    end

    % re-orders points such that point furthest in the x-axis is first on
    % the list of points    
    indx_ = find(points_um(:,1) == max(points_um(:,1)));
    points_um = [points_um(indx_:end,:);points_um(1:indx_-1,:)];

    % finds arc-length of distance around nuclei shape
    dt = diff(points_um);
    dis = zeros(length(dt),1);
    for k = 1:length(dt)
        if k == 1
            dis(k) = norm(dt(k,:));
        else
            dis(k) = dis(k-1) + norm(dt(k,:));
        end
    end
    
    % prepares arc-length data for Fourier series fit
    dis = [0;dis]; % first point is zero arc-length
    dis = normalize(dis,'range')*2*pi; % normalize from arclength to 2pi
    
    % concatenates x and y points separately 
    points_um_x_NUC = points_um(:,1); % concatenates x
    points_um_y_NUC = points_um(:,2); % concatenates y
    
    % computes Fourier series coefficients
    [a] = fourier_series_fit(dis,points_um_x_NUC,total_number_of_series_nuc);
    [b] = fourier_series_fit(dis,points_um_y_NUC,total_number_of_series_nuc);

    % saves Fourier series coefficients
    save(strcat('a_coeff_NUC_',num2str(cell_num),'.mat'),'a')
    save(strcat('b_coeff_NUC_',num2str(cell_num),'.mat'),'b')

    % Fourier series evaluated at these theta values
    thetas = linspace(0,2*pi,10000);
    
    % Evaluation of Fourier series
    [x_predicted_NUC] = fourier_series_evaluate(a,dis);
    [y_predicted_NUC] = fourier_series_evaluate(b,dis);

    %% check Fourier fits
    
    % % uncomment to see. This will generate a figure for every cell. 
    % figure
    % subplot(1,2,1) % compares cell body
    % plot(points_um_x_CELL,points_um_y_CELL); hold on;
    % scatter(x_predicted_CELL,y_predicted_CELL); hold off;
    % xlabel('\mum')
    % ylabel('\mum')
    % title('cell')
    % axis equal
    % subplot(1,2,2) % compares nuclei
    % plot(points_um_x_NUC,points_um_y_NUC); hold on;
    % scatter(x_predicted_NUC,y_predicted_NUC); hold off;
    %  xlabel('\mum')
    % ylabel('\mum')
    % title('nuclei')
    % axis equal

    cell_num = cell_num + 1;

end

%% saves morphological data
SAVEDATA = table((1:cell_num-1)',cell_area,cell_roundness,cell_width,cell_length,cell_aspect_ratio,cell_solidity,...
DAPI_area,DAPI_roundness,DAPI_width,DAPI_length,DAPI_aspect_ratio,DAPI_solidity,activated,DOA);

writetable(SAVEDATA,"DATA.xlsx",'WriteVariableNames',true) 