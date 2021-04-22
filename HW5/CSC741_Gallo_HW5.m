% Created by: Jessica Gallo
% Date Created: 3/8/2021
% Last Modified: 3/8/2021
% CSC741 Digital Image Processing
% HW 5

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

% =========================================================================
% Part 1:
% Generate and add 4 different noise (including Gaussian salt-and-pepper
% (impulse)) noise to 3 images. You must be able to specify the noise mean
% and variance (for example=10, 0.8)
% =========================================================================

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

% =========================================================================
% Part 2
% Estimate noise parameters from a single image
% =========================================================================

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

% =========================================================================
% Part 3
% Apply a class of Mean (Arithmetic, weighted mean, Geometric, Harmonic,
% mean filter)
% =========================================================================
% Arithmetic Mean Filter
A = fspecial('average', [7,7]);
Im1G_denoise1_AM = filter2(A,Im1G_noise1);
Im1G_denoise2_AM = filter2(A,Im1G_noise2);
Im1G_denoise3_AM = filter2(A,Im1G_noise3);
Im1G_denoise4_AM = filter2(A,Im1G_noise4);

figure(11)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_AM)); xlabel('Denoise Gaussian'); title('Arithmetic Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_AM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_AM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_AM)); xlabel('Denoise Speckle');

figure(15)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_AM)); ylabel('Denoise Gaussian'); title('Arithmetic Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_AM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_AM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_AM)); ylabel('Denoise Speckle');

arithmean1 = geomean(Im1G_denoise1_AM);
fprintf('\nThe arithmetic mean of gaussian image is: %0.4d', arithmean1);
arithmean2 = geomean(Im1G_denoise2_AM);
fprintf('\nThe arithmetic mean of salt&pepper image is: %0.4d', arithmean2);
arithmean3 = geomean(Im1G_denoise3_AM);
fprintf('\nThe arithmetic mean of poisson image is: %0.4d', arithmean3);
arithmean4 = geomean(Im1G_denoise4_AM);
fprintf('\nThe arithmetic mean of speckle image is: %0.4d', arithmean4);

% Weighted Mean Filter
f2=1/16*[1,2,1;2,4,2;1,2,1];
Im1G_denoise1_WM = filter2(f2, Im1G_noise1); % filter2 function
Im1G_denoise2_WM = filter2(f2, Im1G_noise2);
Im1G_denoise3_WM = filter2(f2, Im1G_noise3);
Im1G_denoise4_WM = filter2(f2, Im1G_noise4);

figure(12)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(uint8(Im1G_denoise1_WM)); xlabel('Denoise Gaussian'); title('Weighted Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_WM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_WM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_WM)); xlabel('Denoise Speckle');

figure(16)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussia4n'); title('Original Noise');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_WM)); ylabel('Denoise Gaussian'); title('Weighted Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_WM)); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_WM)); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_WM)); ylabel('Denoise Speckle');

% Geometric Mean Filter
h = ones(5,5);
Im1G_denoise1_GM = imfilter(log(double(Im1G_noise1)), h, 'replicate');
Im1G_denoise1_GM = exp(Im1G_denoise1_GM);
Im1G_denoise1_GM = Im1G_denoise1_GM .^ (1/numel(h));
Im1G_denoise2_GM = imfilter(log(double(Im1G_noise2)), h, 'replicate');
Im1G_denoise2_GM = exp(Im1G_denoise2_GM);
Im1G_denoise2_GM = Im1G_denoise2_GM .^ (1/numel(h));
Im1G_denoise3_GM = imfilter(log(double(Im1G_noise3)), h, 'replicate');
Im1G_denoise3_GM = exp(Im1G_denoise3_GM);
Im1G_denoise3_GM = Im1G_denoise3_GM .^ (1/numel(h));
Im1G_denoise4_GM = imfilter(log(double(Im1G_noise4)), h, 'replicate');
Im1G_denoise4_GM = exp(Im1G_denoise4_GM);
Im1G_denoise4_GM = Im1G_denoise4_GM .^ (1/numel(h));

figure(13)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imshow(Im1G_denoise1_GM, []); xlabel('Denoise Gaussian'); title('Geomtric Mean');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(Im1G_denoise2_GM, []); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(Im1G_denoise3_GM, []); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(Im1G_denoise4_GM, []); xlabel('Denoise Speckle');

figure(17)
subplot(4, 2, 1); imhist(Im1G_noise1); ylabel('Gaussian'); title('Original Noise');
subplot(4, 2, 2); imhist(Im1G_denoise1_GM); ylabel('Denoise Gaussian'); title('Geomtric Mean');
subplot(4, 2, 3); imhist(Im1G_noise2); ylabel('Salt&Pepper');
subplot(4, 2, 4); imhist(Im1G_denoise2_GM); ylabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); ylabel('Poisson');
subplot(4, 2, 6); imhist(Im1G_denoise3_GM); ylabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); ylabel('Speckle');
subplot(4, 2, 8); imhist(Im1G_denoise4_GM); ylabel('Denoise Speckle');

denoise1_GM = imhist(Im1G_denoise1_GM);
arithmean1 = geomean(denoise1_GM);
fprintf('\nThe geometric mean of gaussian image is: %0.4d', arithmean1);
denoise2_GM = imhist(Im1G_denoise2_GM);
arithmean2 = geomean(denoise2_GM);
fprintf('\nThe geometric mean of salt&pepper image is: %0.4d', arithmean2);
denoise3_GM = imhist(Im1G_denoise3_GM);
geomean3 = geomean(denoise3_GM);
fprintf('\nThe geometric mean of poisson image is: %0.4d', geomean3);
denoise4_GM = imhist(Im1G_denoise4_GM);
arithmean4 = geomean(denoise4_GM);
fprintf('\nThe geometric mean of speckle image is: %0.4d', arithmean4);

% Harmonic Mean Filter
Im1G_noise1D = double(Im1G_noise1);
Im1G_noise2D = double(Im1G_noise2);
Im1G_noise3D = double(Im1G_noise3);
Im1G_noise4D = double(Im1G_noise4);
Im1G_denoise1_HM = (9*9)./imfilter(1./(Im1G_noise1D+eps), ones(3, 3), 'replicate');
Im1G_denoise2_HM = (9*9)./imfilter(1./(Im1G_noise2D+eps), ones(3, 3), 'replicate');
Im1G_denoise3_HM = (9*9)./imfilter(1./(Im1G_noise3D+eps), ones(3, 3), 'replicate');
Im1G_denoise4_HM = (9*9)./imfilter(1./(Im1G_noise4D+eps), ones(3, 3), 'replicate');

figure(14)
subplot(4, 2, 1); imshow(Im1G_noise1); xlabel('Gaussian'); title('Harmonic Mean');
subplot(4, 2, 2); imshow(uint16(Im1G_denoise1_HM)); xlabel('Denoise Gaussian');
subplot(4, 2, 3); imshow(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imshow(uint8(Im1G_denoise2_HM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imshow(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imshow(uint8(Im1G_denoise3_HM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imshow(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imshow(uint8(Im1G_denoise4_HM)); xlabel('Denoise Speckle');

figure(18)
subplot(4, 2, 1); imhist(Im1G_noise1); xlabel('Gaussian'); title('Harmonic Mean');
subplot(4, 2, 2); imhist(uint8(Im1G_denoise1_HM)); xlabel('Denoise Gaussian');
subplot(4, 2, 3); imhist(Im1G_noise2); xlabel('Salt&Pepper');
subplot(4, 2, 4); imhist(uint8(Im1G_denoise2_HM)); xlabel('Denoise Salt&Pepper');
subplot(4, 2, 5); imhist(Im1G_noise3); xlabel('Poisson');
subplot(4, 2, 6); imhist(uint8(Im1G_denoise3_HM)); xlabel('Denoise Poisson');
subplot(4, 2, 7); imhist(Im1G_noise4); xlabel('Speckle');
subplot(4, 2, 8); imhist(uint8(Im1G_denoise4_HM)); xlabel('Denoise Speckle');


% Part 4
% =========================================================================
%g(x,y) = input image
%f(x,y)=degraded image
%H: degradation process
%n(x,y) = noise



