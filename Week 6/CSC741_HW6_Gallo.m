% Jessica Gallo
% CSC741 Digital Image Processing
% Prof. Sos Agaian
% Lecture 6 HW

close all; clear;

% Load images
Im1 = imread('thermalImage3.jpg');
Im2 = imread('thermalImage4.jpg');
Im3 = imread('thermalImage5.jpg');
% RBG to Grey
Im1G = rgb2gray(Im1);
Im2G = rgb2gray(Im2);
Im3G = rgb2gray(Im3);
% Resize
Im1G = imresize(Im1G, [120, 120]);
Im2G = imresize(Im2G, [120, 120]);
Im3G = imresize(Im3G, [120, 120]);
figure(1)
subplot(3, 1, 1); imshow(Im1G); xlabel('Image 1'); title('Images Greyscale and Resized');
subplot(3, 1, 2); imshow(Im2G); xlabel('Image 2');
subplot(3, 1, 3); imshow(Im3G); xlabel('Image 3');

% 1
% Generate and add 4 different noise (including Gaussian, Salta-and-Pepper (impulse)) noise to 3 images. 
% You must be able to specify the noise mean and variance (for example=10, 0.8)
% ================================================================================================================
% Image 1
Im1G_noise1 = imnoise(Im1G, 'gaussian', 0.02);
Im1G_noise2 = imnoise(Im1G, 'salt & pepper', 0.02); % 0.5 is default
Im1G_noise3 = imnoise(Im1G, 'poisson');
Im1G_noise4 = imnoise(Im1G, 'speckle', 0.02); % 0.5 default

figure(2)
subplot(3, 2, 1); imshow(Im1); xlabel('Original'); title('Noise in Image 1');
subplot(3, 2, 2); imshow(Im1G_noise1); xlabel('Gaussian Noise');
subplot(3, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(3, 2, 4); imshow(Im1G_noise3); xlabel('Poisson');
subplot(3, 2, 5); imshow(Im1G_noise4); xlabel('Speckle');

figure(3)
subplot(3, 2, 1); imhist(Im1); ylabel('Original'); title('Noise in Image 1');
subplot(3, 2, 2); imhist(Im1G_noise1); ylabel('Gaussian Noise');
subplot(3, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(3, 2, 4); imhist(Im1G_noise3); ylabel('Poisson');
subplot(3, 2, 5); imhist(Im1G_noise4); ylabel('Speckle');

% Image2
Im2G_noise1 = imnoise(Im2G, 'gaussian', 0.02);
Im2G_noise2 = imnoise(Im2G, 'salt & pepper', 0.02); % 0.5 is default
Im2G_noise3 = imnoise(Im2G, 'poisson');
Im2G_noise4 = imnoise(Im2G, 'speckle', 0.02); % 0.5 default

figure(4)
subplot(3, 2, 1); imshow(Im2G); xlabel('Original'); title('Noise in Image 2');
subplot(3, 2, 2); imshow(Im2G_noise1); xlabel('Gaussian Noise');
subplot(3, 2, 3); imshow(Im2G_noise2); xlabel('Salt&Pepper');
subplot(3, 2, 4); imshow(Im2G_noise3); xlabel('Poisson');
subplot(3, 2, 5); imshow(Im2G_noise4); xlabel('Speckle');

figure(5)
subplot(3, 2, 1); imhist(Im2G); ylabel('Original'); title('Noise in Image 2');
subplot(3, 2, 2); imhist(Im2G_noise1); ylabel('Gaussian Noise');
subplot(3, 2, 3); imhist(Im2G_noise2); ylabel('Salt&Pepper');
subplot(3, 2, 4); imhist(Im2G_noise3); ylabel('Poisson');
subplot(3, 2, 5); imhist(Im2G_noise4); ylabel('Speckle');

% Image 3
Im3G_noise1 = imnoise(Im3G, 'gaussian', 0.02); % variance = 0.4 | mean = zeromean
Im3G_noise2 = imnoise(Im3G, 'salt & pepper', 0.02); % variance = 0.5 is default
Im3G_noise3 = imnoise(Im3G, 'poisson');
Im3G_noise4 = imnoise(Im3G, 'speckle', 0.02); % variance = 0.5 default

figure(6)
subplot(3, 2, 1); imshow(Im3G); xlabel('Original'); title('Noise in Image 3');
subplot(3, 2, 2); imshow(Im3G_noise1); xlabel('Gaussian Noise');
subplot(3, 2, 3); imshow(Im3G_noise2); xlabel('Salt&Pepper');
subplot(3, 2, 4); imshow(Im3G_noise3); xlabel('Poisson');
subplot(3, 2, 5); imshow(Im3G_noise4); xlabel('Speckle');

figure(7)
subplot(3, 2, 1); imhist(Im3G); ylabel('Original'); title('Noise in Image 3');
subplot(3, 2, 2); imhist(Im3G_noise1); ylabel('Gaussian Noise');
subplot(3, 2, 3); imhist(Im3G_noise2); ylabel('Salt&Pepper');
subplot(3, 2, 4); imhist(Im3G_noise3); ylabel('Poisson');
subplot(3, 2, 5); imhist(Im3G_noise4); ylabel('Speckle');

% 2
% Estimate noise parameters from a single image
% ================================================================================================================
% Segment a small patch of image with constant grey level
Im1G_noise1_GC = grayconnected(Im1G_noise1,50,50); % segments the image at coordinates 50,50
figure(8)
subplot(2, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian Noise Image');
subplot(2, 2, 2); imhist(Im1G_noise1); ylabel('Gassian Noise Histogram');
subplot(2, 2, 3); imshow(labeloverlay(Im1G_noise1,Im1G_noise1_GC)); xlabel('Overlapped Image'); % shows the original image with the new image
subplot(2, 2, 4); imhist(labeloverlay(Im1G_noise1,Im1G_noise1_GC)); ylabel('Overlapped Histogram'); % histogram of overlapped image
% Inspect histogram
figure(9)
subplot(3, 1, 1); imshow(Im1G); xlabel('Original');
subplot(3, 1, 2); imshow(Im1G_noise1); xlabel('Image 1 Noise (Gaussian)');
subplot(3, 1, 3); imshow(Im1G_noise1_GC); xlabel('Small Patch of Noisy Image (Gaussian)');
% Compute the mean and variance
x = imhist(Im1G_noise1_GC);
Im1G_noise1_mean = mean2(x);
fprintf('The mean of the small patch of the noise image is %0.4d', Im1G_noise1_mean);
Im1G_noise1_var = var(x);
fprintf('\nThe variance of the small patch of the noise image is %0.4d', Im1G_noise1_var);
Im1G_noise1_std = std(x);
fprintf('\nThe standard deviation of the small patch of the noise image is %0.4d', Im1G_noise1_std);

% 3
% Apply a class of order-statistic (median, max/min, weighted median, midpoint, alpha-trimmed mean filter) filters
% on generated noisy images
% ================================================================================================================

% Median
Im1G_denoise1_M = medfilt2(Im1G_noise1);
Im1G_denoise2_M = medfilt2(Im1G_noise2);
Im1G_denoise3_M = medfilt2(Im1G_noise3);
Im1G_denoise4_M = medfilt2(Im1G_noise4);

figure(10)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_M)); xlabel('Denoise Gaussian'); title('Median Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_M)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_M)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_M)); xlabel('Denoise Speckle');

figure(11)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_M)); ylabel('Denoise Gaussian'); title('Median Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_M)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_M)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_M)); ylabel('Denoise Speckle');

% Max/Min
Im1G_denoise1_MaxM = imdilate(Im1G_noise1);
Im1G_denoise2_MaxM = imdilate(Im1G_noise2);
Im1G_denoise3_MaxM = imdilate(Im1G_noise3);
Im1G_denoise4_MaxM = imdilate(Im1G_noise4);

Im1G_denoise1_MinM = imerode(Im1G_noise1);
Im1G_denoise2_MinM = imerode(Im1G_noise2);
Im1G_denoise3_MinM = imerode(Im1G_noise3);
Im1G_denoise4_MinM = imerode(Im1G_noise4);

figure(12)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_MaxM)); xlabel('Denoise Gaussian'); title('Max Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_MaxM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_MaxM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_MaxM)); xlabel('Denoise Speckle');
figure(13)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_MinM)); xlabel('Denoise Gaussian'); title('Min Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_MinM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_MinM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_MinM)); xlabel('Denoise Speckle');

figure(14)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_MaxM)); ylabel('Denoise Gaussian'); title('Max Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_MaxM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_MaxM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_MaxM)); ylabel('Denoise Speckle');
figure(15)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_MinM)); ylabel('Denoise Gaussian'); title('Min Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_MinM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_MinM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_MinM)); ylabel('Denoise Speckle');
%{
% Weighted Median
figure(14)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_AM)); xlabel('Denoise Gaussian'); title('Weighted Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_AM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_AM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_AM)); xlabel('Denoise Speckle');

figure(15)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_AM)); ylabel('Denoise Gaussian'); title('Weighted Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_AM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_AM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_AM)); ylabel('Denoise Speckle');

% Midpoint
figure(16)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_AM)); xlabel('Denoise Gaussian'); title('Midpoint Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_AM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_AM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_AM)); xlabel('Denoise Speckle');

figure(17)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_AM)); ylabel('Denoise Gaussian'); title('Midpoint Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_AM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_AM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_AM)); ylabel('Denoise Speckle');

% Alpha-Trimmed Mean
AT_mask = [1 1 1; 1 1 1; 1 1 1];


figure(18)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_AM)); xlabel('Denoise Gaussian'); title('Alpha-Trimmed Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_AM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_AM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_AM)); xlabel('Denoise Speckle');

figure(19)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_AM)); ylabel('Denoise Gaussian'); title('Alpha-Trimmed Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_AM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_AM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_AM)); ylabel('Denoise Speckle');
%}