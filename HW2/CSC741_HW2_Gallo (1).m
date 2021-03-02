close all; clear;

% Image Preprocessing =====================================================
% Defining Images
Im1 = imread('thermalImage3.jpg');
%Im2 = imread('thermalImage2.jpg');

% RGB to Grey
Im1G = rgb2gray(Im1); 
%Im2G = rgb2gray(Im2);

figure(1)
subplot(2, 1, 1); imshow(Im1); xlabel('Image 1 Original');
%subplot(2, 2, 2); imshow(Im2); xlabel('Image 2 Original');
subplot(2, 1, 2); imshow(Im1G); xlabel('Image 1 Grey');
%subplot(2, 2, 4); imshow(Im2G); xlabel('Image 2 Grey');

figure(1)
% Resize to same size
Im1G = imresize(Im1G, [120, 120]);
%Im2G = imresize(Im2G, [500, 700]);

subplot(2, 1, 1); imshow(Im1G); xlabel('Image 1 Resized');
%subplot(2, 1, 2); imshow(Im2G); xlabel('Image 2 Resized');

% Image subtraction from last week
%subIm = imsubtract(Im1G,Im2G);

%figure(12);
%subplot(1, 1, 1); imshow(subIm); xlabel('Image Subtraction');

% Histograms ==============================================================
% Bar histogram with greyscale image vectorized
figure(2)
Im1G_vector = Im1G(:); % converts to a single vector
subplot(2, 1, 1); imhist(Im1G); ylabel('Grey Image Histogram');
%subplot(2, 1, 2); imhist(Im1G_vector); ylabel('Grey Image 1 Vector Histogram');

% Histogram sliding % Contrast Stretching
% ----
% adding/subtracting a constant brightness value to all pixels in the image
% ----
Im1G_HSlide=im2double(Im1G);
bright_add = 0.2;
h2=Im1G_HSlide+bright_add;

%B = 2;
%h3 = Im1G_vector;
%Im1G_contrast = (2^B - 1)*(Im1G - min(Im1G_vector)) / (max(Im1G_vector)-min(Im1G_vector));
Im1G_contrast = imadjust(Im1G, stretchlim(Im1G, [0.05, 0.95]),[]);

figure(3);
subplot(3, 2, 1); imshow(Im1G); xlabel('Original Image');
subplot(3, 2, 2); imhist(Im1G); ylabel('Original Histogram');
subplot(3, 2, 3); imshow(h2), xlabel('Histogram Slide Image');
subplot(3, 2, 4); imhist(h2), ylabel('Histogram Slide Histogram');
subplot(3, 2, 5); imshow(uint8(Im1G_contrast)); xlabel('Contrast Stretching')
subplot(3, 2, 6); imhist(Im1G_contrast); ylabel('Contrast Stretching');

% Contrast Enhancement ====================================================
% ----
% An image enhancement technique that improves the contrast in an image by
% expanding the dynamic range of intensity values it contains
% ----

% Adjust & Adapted Stretching & Histogram Stretching
Im1G_adjust = imadjust(Im1G);
Im1G_HStretch = histeq(Im1G);
Im1G_adapthisteq = adapthisteq(Im1G);

figure(4)
title("Original and imadjust, histeq, and adapthisteq")
subplot(2, 2, 1); imshow(Im1G); ylabel('Original');
subplot(2, 2, 2); imshow(Im1G_adjust); ylabel('Adjusted');
subplot(2, 2, 3); imshow(Im1G_HStretch); ylabel('Histogram Stretch');
subplot(2, 2, 4); imshow(Im1G_adapthisteq); ylabel('Adapted Stretch');

figure(5)
subplot(2, 2, 1); imhist(Im1G); ylabel('Original');
subplot(2, 2, 2); imhist(Im1G_adjust); ylabel('Adjusted');
subplot(2, 2, 3); imhist(Im1G_HStretch); ylabel('Histogram Stretch');
subplot(2, 2, 4); imhist(Im1G_adapthisteq); ylabel('Adapted Stretch');

% Noisy Image
ImG_noisy = imnoise(Im1G, 'gaussian', 0.04);

figure(6)
subplot(2, 2, 1); imshow(Im1G); xlabel('Original');
subplot(2, 2, 2); imhist(Im1G); ylabel('Original Histogram');
subplot(2, 2, 3); imshow(ImG_noisy); xlabel('Noisy Image');
subplot(2, 2, 4); imhist(ImG_noisy); ylabel('Noisy Histogram');

% MSE (Mean Squared Error)
mse1 = immse(ImG_noisy, Im1G);
fprintf('The MSE between the original image and the noise image is %0.4f', mse1);
Im1G_double2 = double(Im1G);
mse2 = immse(h2, Im1G_double2);
fprintf('\nThe MSE between the original image and the image with histogram sliding is %0.4f', mse2);
mse3 = immse(Im1G_contrast, Im1G);
fprintf('\nThe MSE between the original image and the image with constrast stretching is %0.4f', mse3);
mse4 = immse(Im1G_adjust, Im1G);
fprintf('\nThe MSE between the original image and the image with adjusted contrast is %0.4f', mse4);
mse5 = immse(Im1G_HStretch, Im1G);
fprintf('\nThe MSE between the original image and the image with histogram stretch is %0.4f', mse5);
mse6 = immse(Im1G_adapthisteq, Im1G);
fprintf('\nThe MSE between the original image and the image with adapted histogram stretch is %0.4f', mse6);

fprintf('\n-----------------------------------------------------------------------------------------------');

% PSNR (Peak Signal to Noise Ratio)
[peaksnr1, snr1] = psnr(ImG_noisy, Im1G);
fprintf('\nThe PSNR between the original image and the noisy image is %0.4f', peaksnr1);
fprintf('\nThe SNR between the original image and the noisy image is %0.4f', snr1);
[peaksnr2, snr2] = psnr(h2, Im1G_double2);
fprintf('\nThe PSNR between the original image and the image with histogram sliding is %0.4f', peaksnr2);
fprintf('\nThe SNR between the original image and the image with histogram sliding is %0.4f', snr2);
[peaksnr3, snr3] = psnr(Im1G_contrast, Im1G);
fprintf('\nThe PSNR between the original image and the image with constrast stretching is %0.4f', peaksnr3);
fprintf('\nThe SNR between the original image and the image constrast stretching is %0.4f', snr3);
[peaksnr4, snr4] = psnr(Im1G_adjust, Im1G);
fprintf('\nThe PSNR between the original image and the image with adjusted contrast is %0.4f', peaksnr4);
fprintf('\nThe SNR between the original image and the image with adjusted contrast is %0.4f', snr4);
[peaksnr5, snr5] = psnr(Im1G_HStretch, Im1G);
fprintf('\nThe PSNR between the original image and the image with histogram stretch is %0.4f', peaksnr5);
fprintf('\nThe SNR between the original image and the image with histogram stretch is %0.4f', snr5);
[peaksnr6, snr6] = psnr(Im1G_adapthisteq, Im1G);
fprintf('\nThe PSNR between the original image and the image adapted histogram is %0.4f', peaksnr6);
fprintf('\nThe SNR between the original image and the image adapted histogram is %0.4f', snr6);

fprintf('\n-----------------------------------------------------------------------------------------------');

% Euclidean Distance
% Original Grey Image 1 vs Noisy Image
Im1G_noisy_vector = ImG_noisy(:);
edistance1 = sqrt(sum((Im1G_vector-Im1G_noisy_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the Noisy Image is %0.4f', edistance1);
Im1G_vector2 = Im1G_double2(:);
Im1G_h2_vector = h2(:);
edistance2 = sqrt(sum((Im1G_vector2-Im1G_h2_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the histogram sliding Image is %0.4f', edistance2);
Im1G_contrast_vector = Im1G_contrast(:);
edistance3 = sqrt(sum((Im1G_vector-Im1G_contrast_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the constrast stretching Image is %0.4f', edistance3);
Im1G_adjust_vector = Im1G_adjust(:);
edistance4 = sqrt(sum((Im1G_vector-Im1G_adjust_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the adjusted contrast Image is %0.4f', edistance4);
Im1G_HStretch_vector = Im1G_HStretch(:);
edistance5 = sqrt(sum((Im1G_vector-Im1G_HStretch_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the histogram stretch Image is %0.4f', edistance5);
Im1G_adapthisteq_vector = Im1G_adapthisteq(:);
edistance6 = sqrt(sum((Im1G_vector-Im1G_adapthisteq_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the adapted histogram Image is %0.4f', edistance6);

% Manhattan Distance
% Original Grey Image
ImG_CB1 = bwdist(Im1G, 'cityblock');
ImG_CB2 = bwdist(ImG_noisy, 'cityblock');
ImG_CB3 = bwdist(h2, 'cityblock');
ImG_CB4 = bwdist(Im1G_contrast, 'cityblock');
ImG_CB5 = bwdist(Im1G_adjust, 'cityblock');
ImG_CB6 = bwdist(Im1G_HStretch, 'cityblock');
ImG_CB7 = bwdist(Im1G_adapthisteq, 'cityblock');
figure(7)
subplot(2, 2, 1); imshow(ImG_CB1); xlabel('Original');
subplot(2, 2, 2); imshow(ImG_CB2); xlabel('Noisy');
subplot(2, 2, 3); imshow(ImG_CB3); xlabel('Histogram Sliding');
subplot(2, 2, 4); imshow(ImG_CB4); xlabel('Contrast Stretching');
figure(10);
subplot(2, 2, 1); imshow(ImG_CB5); xlabel('Adjust Contrast');
subplot(2, 2, 2); imshow(ImG_CB6); xlabel('Histogram Stretching');
subplot(2, 2, 3); imshow(ImG_CB7); xlabel('Adapted Histogram Streching');

% Gamma Correction
gamma = 1.05;
Im1G_double = double(Im1G);
g1_gamma = Im1G_double.^gamma;

figure(8);
subplot(2, 1, 1); imshow(Im1G); xlabel('Original');
subplot(2, 1, 2); imshow(uint8(g1_gamma)); xlabel('Gamma Correction');