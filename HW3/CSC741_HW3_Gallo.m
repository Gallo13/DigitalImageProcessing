% Jessica Gallo
% Created: 3/1/2021
% Last Modified: 3/4/2021
% CSC741 Digital Image Processing
% Week 4 HW

close all; clear;

% Image Preprocessing =====================================================
% load images
Im1 = imread('thermalImage3.jpg');
Im2 = imread('thermalImage4.jpg');
% images to greyscale
Im1G = rgb2gray(Im1);
Im2G = rgb2gray(Im2);
% images resized to same size
Im1G = imresize(Im1G, [120, 120]);
Im2G = imresize(Im2G, [120, 120]);
figure(1)
subplot(2, 1, 1); imshow(Im1G); title('Images Greyscaled and Resized');
subplot(2, 1, 2); imshow(Im2G);

% Gamma Correction ========================================================
Im1G_double = double(Im1G);
Im2G_double = double(Im2G);

gamma1 = 1.05;
g1_gamma1 = Im1G_double.^gamma1;
g2_gamma1 = Im2G_double.^gamma1;

gamma2 = 0.95;
g1_gamma2 = Im1G_double.^gamma2;
g2_gamma2 = Im2G_double.^gamma2;

figure(2);
subplot(3, 2, 1); imshow(Im1G); xlabel('Original 1'); title('Gamma Correction Image 1');
subplot(3, 2, 3); imshow(uint8(g1_gamma1)); xlabel('Gamma Correction > 1');
subplot(3, 2, 5); imshow(uint8(g1_gamma2)); xlabel('Gamma Correction < 1');
subplot(3, 2, 2); imshow(Im2G); xlabel('Original 2'); title('Gamma Correction Image 2');
subplot(3, 2, 4); imshow(uint8(g2_gamma1)); xlabel('Gamma Correction > 1');
subplot(3, 2, 6); imshow(uint8(g2_gamma2)); xlabel('Gamma Correction < 1');

% Alpha Blending ==========================================================
%Alpha blending is a convex combination of two colors allowing for 
%transparency effects in computer graphics. The value of alpha in the color 
%code ranges from 0.0 to 1.0, where 0.0 represents a fully transparent 
%color, and 1.0 represents a fully opaque color.

%The value of the resulting color when color Value1 with an alpha value of 
%Alpha is drawn over a background of color Value0 is given by:

%Value = Value0(1.0 - Alpha) + Value1(Alpha)

alpha = 0.6;
alphaBlendIm = alpha * Im1G + (1 - alpha) * Im2G;
imwrite(alphaBlendIm, 'alphaBlendIm.jpg');

figure(3)
subplot(3, 1, 1); imshow(Im1G); xlabel('Image 1'); title('Alpha Blending');
subplot(3, 1, 2); imshow(Im2G); xlabel('Image 2');
subplot(3, 1, 3); imshow(alphaBlendIm); xlabel('Alpha Blended Image');

% Histograms ==============================================================
% Bar histogram with greyscale image vectorized
figure(4)
Im1G_vector = Im1G(:); % converts to a single vector
Im2G_vector = Im2G(:);
subplot(2, 1, 1); imhist(Im1G_vector); ylabel('Grey Image 1 Histogram');
subplot(2, 1, 2); imhist(Im2G_vector); ylabel('Grey Image 2 Histogram');

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

figure(5);
subplot(3, 2, 1); imshow(Im1G); xlabel('Original Image');
subplot(3, 2, 2); imhist(Im1G); ylabel('Original Histogram');
subplot(3, 2, 3); imshow(h2), xlabel('Histogram Slide Image');
subplot(3, 2, 4); imhist(h2), ylabel('Histogram Slide Histogram');
subplot(3, 2, 5); imshow(uint8(Im1G_contrast)); xlabel('Contrast Stretching')
subplot(3, 2, 6); imhist(Im1G_contrast); ylabel('Contrast Stretching');

% Contrast Enhancement ====================================================

% Adjust & Adapted Stretching & Histogram Stretching
Im1G_adjust = imadjust(Im1G);
Im1G_HEqualization = histeq(Im1G);
Im1G_adapthisteq = adapthisteq(Im1G);

figure(6)
title("Original and imadjust, histeq, and adapthisteq")
subplot(2, 2, 1); imshow(Im1G); ylabel('Original');
subplot(2, 2, 2); imshow(Im1G_adjust); ylabel('Adjusted');
subplot(2, 2, 3); imshow(Im1G_HEqualization); ylabel('Histogram Equalization');
subplot(2, 2, 4); imshow(Im1G_adapthisteq); ylabel('Adapted Stretch');

figure(7)
subplot(2, 2, 1); imhist(Im1G); ylabel('Original');
subplot(2, 2, 2); imhist(Im1G_adjust); ylabel('Adjusted');
subplot(2, 2, 3); imhist(Im1G_HEqualization); ylabel('Histogram Equalization');
subplot(2, 2, 4); imhist(Im1G_adapthisteq); ylabel('Adapted Stretch');

% Noisy Image
ImG_noisy = imnoise(Im1G, 'gaussian', 0.04);

figure(8)
subplot(2, 2, 1); imshow(Im1G); xlabel('Original');
subplot(2, 2, 2); imhist(Im1G); ylabel('Original Histogram');
subplot(2, 2, 3); imshow(ImG_noisy); xlabel('Noisy Image');
subplot(2, 2, 4); imhist(ImG_noisy); ylabel('Noisy Histogram');

% Distance Measurements ===================================================
g1_gamma1 = uint8(g1_gamma1);
g1_gamma2 = uint8(g1_gamma2);
% MSE (Mean Squared Error)
mse0 = immse(g1_gamma1, Im1G);
fprintf('The MSE between the original image and the gamma > 1 image is %0.4f', mse0);
mse01 = immse(g1_gamma2, Im1G);
fprintf('\nThe MSE between the original image and the gamma < 1 image is %0.4f', mse01);
mse1 = immse(ImG_noisy, Im1G);
fprintf('\nThe MSE between the original image and the noise image is %0.4f', mse1);
Im1G_double2 = double(Im1G);
mse2 = immse(h2, Im1G_double2);
fprintf('\nThe MSE between the original image and the image with histogram sliding is %0.4f', mse2);
mse3 = immse(Im1G_contrast, Im1G);
fprintf('\nThe MSE between the original image and the image with constrast stretching is %0.4f', mse3);
mse4 = immse(Im1G_adjust, Im1G);
fprintf('\nThe MSE between the original image and the image with adjusted contrast is %0.4f', mse4);
mse5 = immse(Im1G_HEqualization, Im1G);
fprintf('\nThe MSE between the original image and the image with histogram equalization is %0.4f', mse5);
mse6 = immse(Im1G_adapthisteq, Im1G);
fprintf('\nThe MSE between the original image and the image with adapted histogram stretch is %0.4f', mse6);

fprintf('\n-----------------------------------------------------------------------------------------------');

% PSNR (Peak Signal to Noise Ratio)
[peaksnr0, snr0] = psnr(g1_gamma1, Im1G);
fprintf('\nThe PSNR between the original image and the gamma > 1 image is %0.4f', peaksnr0);
fprintf('\nThe SNR between the original image and the gamma > 1 image is %0.4f', snr0);
[peaksnr01, snr01] = psnr(g1_gamma2, Im1G);
fprintf('\nThe PSNR between the original image and the gamma < 1 image is %0.4f', peaksnr01);
fprintf('\nThe SNR between the original image and the gamma < 1 image is %0.4f', snr01);
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
[peaksnr5, snr5] = psnr(Im1G_HEqualization, Im1G);
fprintf('\nThe PSNR between the original image and the image with histogram equalization is %0.4f', peaksnr5);
fprintf('\nThe SNR between the original image and the image with histogram equalization is %0.4f', snr5);
[peaksnr6, snr6] = psnr(Im1G_adapthisteq, Im1G);
fprintf('\nThe PSNR between the original image and the image adapted histogram is %0.4f', peaksnr6);
fprintf('\nThe SNR between the original image and the image adapted histogram is %0.4f', snr6);

fprintf('\n-----------------------------------------------------------------------------------------------');

% Euclidean Distance
g1_gamma1_vector = g1_gamma1(:);
edistance0 = sqrt(sum((Im1G_vector-g1_gamma1_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the Gamma > 1 is %0.4f', edistance0);
g1_gamma2_vector = g1_gamma2(:);
edistance01 = sqrt(sum((Im1G_vector-g1_gamma2_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the Gamma < 1 is %0.4f', edistance01);
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
Im1G_HEqualization_vector = Im1G_HEqualization(:);
edistance5 = sqrt(sum((Im1G_vector-Im1G_HEqualization_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the histogram equalization Image is %0.4f', edistance5);
Im1G_adapthisteq_vector = Im1G_adapthisteq(:);
edistance6 = sqrt(sum((Im1G_vector-Im1G_adapthisteq_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the adapted histogram Image is %0.4f', edistance6);

fprintf('\n-----------------------------------------------------------------------------------------------');

% Manhattan Distance
ImG_CB1 = bwdist(Im1G, 'cityblock');
ImG_CB2 = bwdist(ImG_noisy, 'cityblock');
ImG_CB3 = bwdist(h2, 'cityblock');
ImG_CB4 = bwdist(Im1G_contrast, 'cityblock');
ImG_CB5 = bwdist(Im1G_adjust, 'cityblock');
ImG_CB6 = bwdist(Im1G_HEqualization, 'cityblock');
ImG_CB7 = bwdist(Im1G_adapthisteq, 'cityblock');

%fprintf('\nThe Manhattan Distance of the Original Image is %0.4f', ImG_CB1);
%fprintf('\nThe Manhattan Distance of the noisy image is %0.4f', ImG_CB2);
%fprintf('\nThe Manhattan Distance of the historgeam sliding image is %0.4f', ImG_CB3);
%fprintf('\nThe Manhattan Distance of the contrast stretching image is %0.4f', ImG_CB4);
%fprintf('\nThe Manhattan Distance of the adjusted contrast image is %0.4f', ImG_CB5);
%fprintf('\nThe Manhattan Distance of the histogram stretch image is %0.4f', ImG_CB6);
%fprintf('\nThe Manhattan Distance of the adapted histogram image is %0.4f', ImG_CB7);

figure(9)
subplot(2, 2, 1); imshow(ImG_CB1); xlabel('Original');
subplot(2, 2, 2); imshow(ImG_CB2); xlabel('Noisy');
subplot(2, 2, 3); imshow(ImG_CB3); xlabel('Histogram Sliding');
subplot(2, 2, 4); imshow(ImG_CB4); xlabel('Contrast Stretching');
figure(10);
subplot(2, 2, 1); imshow(ImG_CB5); xlabel('Adjust Contrast');
subplot(2, 2, 2); imshow(ImG_CB6); xlabel('Histogram Equalization');
subplot(2, 2, 3); imshow(ImG_CB7); xlabel('Adapted Histogram Streching');

% Results Together for Comparison =========================================
figure(11)
subplot(4, 4, 1); imshow(Im1G); xlabel('Original');
subplot(4, 4, 2); imshow(uint8(g1_gamma1)); xlabel('Gamma Correction > 1');
subplot(4, 4, 3); imshow(uint8(g1_gamma2)); xlabel('Gamma Correction < 1');
subplot(4, 4, 4); imshow(h2); xlabel('Histogram Slide');
subplot(4, 4, 5); imshow(Im1G_contrast); xlabel('Contrast Stretching');
subplot(4, 4, 6); imshow(Im1G_adjust); xlabel('Contrast Adjust');
subplot(4, 4, 7); imshow(Im1G_HEqualization); xlabel('Histogram Equalization');
subplot(4, 4, 8); imshow(Im1G_adapthisteq); xlabel('Adapted Histogram Eq.');

% Combination of 2 best versions ==========================================
Im1G_adjustAdHist = adapthisteq(Im1G_adjust);
Im1G_adjContrast = imadjust(Im1G_adjust, stretchlim(Im1G_adjust, [0.05, 0.95]),[]);
Im1G_h2contrast = imadjust(h2, stretchlim(h2, [0.05, 0.95]),[]);
figure(12)
subplot(3, 1, 1); imshow(Im1G_adjustAdHist); xlabel('Adapted Hist Eq. on Contrast Adjust Image');
subplot(3, 1, 2); imshow(Im1G_adjContrast); xlabel('Contrast Stretch on Contrast Adjust Image');
subplot(3, 1, 3); imshow(Im1G_h2contrast); xlabel('Contrast Adjust on Histogram Sliding Image');

mse7 = immse(Im1G_adjustAdHist, Im1G);
fprintf('\nThe MSE between the original image and the Adapted Hist Eq. on Contrast Adjust image is %0.4f', mse7);
Im1G_adjContrast = double(Im1G_adjContrast);
mse8 = immse(Im1G_adjContrast, Im1G_double2);
fprintf('\nThe MSE between the original image and the Contrast Stretch on Contrast Adjust image is %0.4f', mse8);
Im1G_h2contrast = double(Im1G_h2contrast);
mse9 = immse(Im1G_h2contrast, Im1G_double2);
fprintf('\nThe MSE between the original image and the Contrast Adjust on Histogram Sliding image is %0.4f', mse9);

fprintf('\n--------------------------')

[peaksnr7, snr7] = psnr(Im1G_adjustAdHist, Im1G);
fprintf('\nThe PSNR between the original image and the Adapted Hist Eq. on Contrast Adjust image is %0.4f', peaksnr7);
fprintf('\nThe SNR between the original image and the Adapted Hist Eq. on Contrast Adjust image is %0.4f', snr7);
[peaksnr8, snr8] = psnr(Im1G_adjContrast, Im1G_double2);
fprintf('\nThe PSNR between the original image and the Contrast Stretch on Contrast Adjust image is %0.4f', peaksnr8);
fprintf('\nThe SNR between the original image and the Contrast Stretch on Contrast Adjust image is %0.4f', snr8);
[peaksnr9, snr9] = psnr(Im1G_h2contrast, Im1G_double2);
fprintf('\nThe PSNR between the original image and Contrast Adjust on Histogram Sliding image is %0.4f', peaksnr9);
fprintf('\nThe SNR between the original image and the Contrast Adjust on Histogram Sliding image is %0.4f', snr9);

fprintf('\n--------------------------')

Im1G_adjustAdHist_vector = Im1G_adjustAdHist(:);
edistance7 = sqrt(sum((Im1G_vector-Im1G_adjustAdHist_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the Adapted Hist Eq. on Contrast Adjust image is %0.4f', edistance7);
Im1G_adjContrast_vector = Im1G_adjContrast(:);
edistance8 = sqrt(sum((Im1G_vector2-Im1G_adjContrast_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the Contrast Stretch on Contrast Adjust image is %0.4f', edistance8);
Im1G_h2contrast_vector = Im1G_h2contrast(:);
edistance9 = sqrt(sum((Im1G_vector2-Im1G_h2contrast_vector).^2));
fprintf('\nThe Euclidean Distance between Original Image 1 and the Contrast Adjust on Histogram Sliding image is %0.4f', edistance9);

