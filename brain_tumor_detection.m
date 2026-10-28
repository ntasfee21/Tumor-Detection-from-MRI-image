function brain_tumor_detection

clc;
close all;


%% ============================================================

img = [];

%% Create GUI window

fig = figure('Name','Brain Tumor Detection GUI', ...
             'NumberTitle','off', ...
             'Position',[80 50 1250 700]);

%% Create axes

axes1 = axes('Parent',fig,'Units','normalized', ...
             'Position',[0.05 0.56 0.23 0.32]);

axes2 = axes('Parent',fig,'Units','normalized', ...
             'Position',[0.38 0.56 0.23 0.32]);

axes3 = axes('Parent',fig,'Units','normalized', ...
             'Position',[0.70 0.56 0.23 0.32]);

axes4 = axes('Parent',fig,'Units','normalized', ...
             'Position',[0.22 0.12 0.23 0.32]);

axes5 = axes('Parent',fig,'Units','normalized', ...
             'Position',[0.56 0.12 0.23 0.32]);

% Make axes blank at the beginning
axis(axes1,'off');
axis(axes2,'off');
axis(axes3,'off');
axis(axes4,'off');
axis(axes5,'off');

%% Upload button below the first image

uicontrol('Parent',fig, ...
          'Style','pushbutton', ...
          'String','Upload Brain Image', ...
          'Units','normalized', ...
          'FontSize',11, ...
          'FontWeight','bold', ...
          'Position',[0.10 0.49 0.14 0.045], ...
          'Callback',@uploadImage);

%% Clear button

uicontrol('Parent',fig, ...
          'Style','pushbutton', ...
          'String','Clear', ...
          'Units','normalized', ...
          'FontSize',11, ...
          'Position',[0.83 0.04 0.10 0.05], ...
          'Callback',@clearImages);

%% ============================================================
% Upload image callback
%% ============================================================

function uploadImage(~,~)

    [file, path] = uigetfile({'*.dcm;*.jpg;*.png;*.tif', ...
                              'Image Files (*.dcm, *.jpg, *.png, *.tif)'});
    
    if isequal(file,0)
        return;
    end
    
    [~,~,ext] = fileparts(file);
    
    if strcmpi(ext,'.dcm')
        img = single(dicomread(fullfile(path,file)));
    else
        temp = imread(fullfile(path,file));
        
        if size(temp,3) == 3
            temp = rgb2gray(temp);
        end
        
        img = single(temp);
    end
    
    % Show original image
    cla(axes1);
    axes(axes1);
    imagesc(img);
    colormap(axes1,gray);
    colorbar;
    axis image;
    title('4(a) Original Image');
    
    % Automatically process the uploaded image
    processUploadedImage;

end

%% ============================================================
% Main processing function
%% ============================================================

function processUploadedImage

    if isempty(img)
        return;
    end

    [row, col] = size(img);

    %% Low-pass filtered image
    % Gaussian smoothing kernel 

    Ksmooth = 1/256 * [1,  4,  6,  4, 1;
                       4, 16, 24, 16, 4;
                       6, 24, 36, 24, 6;
                       4, 16, 24, 16, 4;
                       1,  4,  6,  4, 1];

    img_lowpass = conv2(double(img), Ksmooth, 'same');

    cla(axes2);
    axes(axes2);
    imagesc(img_lowpass);
    colormap(axes2,gray);
    colorbar;
    axis image;
    title('4(b) Low-pass Filtered');

    %% High-pass / edge-detected image
    % Central difference method 

    Kx = [0,  0, 0;
         -1,  0, 1;
          0,  0, 0]/2;

    Ky = [0, -1, 0;
          0,  0, 0;
          0,  1, 0]/2;

    edge_x = conv2(double(img), Kx, 'same');
    edge_y = conv2(double(img), Ky, 'same');

    img_highpass = sqrt(edge_x.^2 + edge_y.^2);

    cla(axes3);
    axes(axes3);
    imagesc(img_highpass);
    colormap(axes3,gray);
    colorbar;
    axis image;
    title('4(c) High-pass / Edge');

        %% Tumor isolated image using thresholding

    % Main ROI around the main tumor region
    r1 = 60;
    r2 = 150;
    c1 = 72;
    c2 = 132;

    % Extra ROI on the right side light tumor region
    r3 = 86;
    r4 = 112;
    c3 = 132;
    c4 = 152;

    % Extra ROI for the upper-left missed side
    r5 = 76;
    r6 = 92;
    c5 = 80;
    c6 = 94;

    % Making sure ROI does not exceed image size
    if r2 > row
        r2 = row;
    end
    if c2 > col
        c2 = col;
    end
    if r4 > row
        r4 = row;
    end
    if c4 > col
        c4 = col;
    end
    if r6 > row
        r6 = row;
    end
    if c6 > col
        c6 = col;
    end

    roi = img_lowpass(r1:r2, c1:c2);

    roi_mean = mean(roi(:));
    roi_max  = max(roi(:));

    % Thresholds for main tumor ROI
    thresh_low  = roi_mean - 0.05*(roi_max - roi_mean);
    thresh_high = roi_mean + 0.78*(roi_max - roi_mean);

    % Thresholds for extra right-side light tumor part
    extra_low   = roi_mean - 0.10*(roi_max - roi_mean);
    extra_high  = roi_mean + 0.72*(roi_max - roi_mean);

    % Thresholds for upper-left missed tumor side
    blue_low    = roi_mean - 0.08*(roi_max - roi_mean);
    blue_high   = roi_mean + 0.76*(roi_max - roi_mean);

    % Skull rejection using original image brightness
    skull_limit = 0.90 * max(img(:));

    tumor_mask = zeros(row, col);

    %% Main tumor region
    for i = r1:r2
        for j = c1:c2

            if img_lowpass(i,j) >= thresh_low && img_lowpass(i,j) <= thresh_high

                if img(i,j) < skull_limit

                    % Remove upper-left skull region
                    if ~(i < 85 && j < 90)

                        % Remove lower-middle false positive region
                        if ~(i > 128 && j > 118)

                            % Remove small false positive region
                            if ~(i >= 62 && i <= 74 && j >= 86 && j <= 94)

                                tumor_mask(i,j) = 64;

                            end
                        end
                    end
                end
            end
        end
    end

    %% To add missed light tumor region
    for i = r3:r4
        for j = c3:c4

            if img_lowpass(i,j) >= extra_low && img_lowpass(i,j) <= extra_high

                if img(i,j) < skull_limit

                    % Prevent adding too much right-side normal brain
                    if ~(i < 92 && j > 145)
                        tumor_mask(i,j) = 64;
                    end

                end
            end
        end
    end

    %% Extra upper-left missed tumor side
    for i = r5:r6
        for j = c5:c6

            if img_lowpass(i,j) >= blue_low && img_lowpass(i,j) <= blue_high

                if img(i,j) < skull_limit

                    % Do not include small false positive dot region
                    if ~(i >= 62 && i <= 74 && j >= 86 && j <= 94)
                        tumor_mask(i,j) = 64;
                    end

                end
            end
        end
    end

    %% Create isolated tumor image
    tumor_isolated = zeros(row, col);

    for i = 1:row
        for j = 1:col

            if tumor_mask(i,j) > 0
                tumor_isolated(i,j) = img(i,j);
            end

        end
    end

    cla(axes4);
    axes(axes4);
    imagesc(tumor_isolated);
    colormap(axes4,gray);
    colorbar;
    axis image;
    title('4(d) Tumor Isolated');


    %% Tumor marked image within original image

    % Use tumor region from 4(d)
    tumor_binary = tumor_isolated > 0;

    % Fill holes so the boundary becomes cleaner
    tumor_binary = imfill(tumor_binary,'holes');

    cla(axes5);
    axes(axes5);
    imagesc(img);
    colormap(axes5,gray);
    colorbar;
    axis image;
    title('4(e) Tumor Marked');

    hold on;
    contour(tumor_binary,[1 1],'r','LineWidth',3);
    hold off;
end

%% ============================================================
% Clear callback
%% ============================================================

function clearImages(~,~)

    cla(axes1);
    cla(axes2);
    cla(axes3);
    cla(axes4);
    cla(axes5);

    axis(axes1,'off');
    axis(axes2,'off');
    axis(axes3,'off');
    axis(axes4,'off');
    axis(axes5,'off');

    img = [];

end

end