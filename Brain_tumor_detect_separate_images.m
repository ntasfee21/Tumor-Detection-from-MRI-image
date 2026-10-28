clc;
clear;
close all;


%% ============================================================


%% Step 1: Load original MRI image

img = single(dicomread('Brain_4.dcm'));

[row, col] = size(img);


%% Step 2:  Original image

figure;
imagesc(img);
colormap(gray);
colorbar;
axis image;
title('Original Brain MRI Image');


%% Step 3: Smoothed / Low-pass filtered image

% Gaussian smoothing kernel
Ksmooth = 1/256 * [1,  4,  6,  4, 1;
                   4, 16, 24, 16, 4;
                   6, 24, 36, 24, 6;
                   4, 16, 24, 16, 4;
                   1,  4,  6,  4, 1];

img_lowpass = conv2(double(img), Ksmooth, 'same');

figure;
imagesc(img_lowpass);
colormap(gray);
colorbar;
axis image;
title('Smoothed / Low-pass Filtered Image');


%% Step 4:  Edge-detected / High-pass filtered image

% Central difference kernels
Kx = [0,  0, 0;
     -1,  0, 1;
      0,  0, 0]/2;

Ky = [0, -1, 0;
      0,  0, 0;
      0,  1, 0]/2;

edge_x = conv2(double(img), Kx, 'same');
edge_y = conv2(double(img), Ky, 'same');

img_highpass = sqrt(edge_x.^2 + edge_y.^2);

figure;
imagesc(img_highpass);
colormap(gray);
colorbar;
axis image;
title('Edge-detected / High-pass Filtered Image');
%% Step 5: Tumor isolated image using thresholding

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

% Extra ROI for the upper-left side
r5 = 76;
r6 = 92;
c5 = 80;
c6 = 94;

roi = img_lowpass(r1:r2, c1:c2);

roi_mean = mean(roi(:));
roi_max  = max(roi(:));

% Thresholds for main tumor ROI
thresh_low  = roi_mean - 0.05*(roi_max - roi_mean);
thresh_high = roi_mean + 0.78*(roi_max - roi_mean);

% Thresholds for extra right-side light tumor part
extra_low   = roi_mean - 0.10*(roi_max - roi_mean);
extra_high  = roi_mean + 0.72*(roi_max - roi_mean);

% Thresholds for the upper-left side
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
                        
                        % Remove the small false positive region
                        if ~(i >= 62 && i <= 74 && j >= 86 && j <= 94)
                            
                            tumor_mask(i,j) = 64;
                            
                        end
                    end
                end
            end
        end
    end
end

%% Extra light tumor region (Added for not missing any tumor region inside the image)
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

%% Extra upper-left side(Added for not missing any tumor region inside the image)
for i = r5:r6
    for j = c5:c6
        
        if img_lowpass(i,j) >= blue_low && img_lowpass(i,j) <= blue_high
            
            if img(i,j) < skull_limit
                
                % Do not include the tiny normal brain region
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

figure;
imagesc(tumor_isolated);
colormap(gray);
colorbar;
axis image;
title( 'Tumor Isolated Image');
%% Step 6: Tumor marked image within original image

tumor_binary = tumor_isolated > 0;

% Fill holes so the boundary becomes cleaner
tumor_binary = imfill(tumor_binary,'holes');

figure;
imagesc(img);
colormap(gray);
colorbar;
axis image;
title('Tumor Marked Image within Original Image');

hold on;
contour(tumor_binary,[1 1],'r','LineWidth',3);
hold off;

%% End

disp('Five separate images have been created successfully.');