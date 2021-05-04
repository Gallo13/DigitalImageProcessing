% Created by: Jessica Gallo
% Date created: 4/23/2021
% Last modified: 4/29/2021
% CSC741 Digital Image Processing
% HW10
% Project 1
% Morphological Image Processing

clear;
clc;
close all;

Im1 = imread('thermalImage3.jpg');
Im1 = rgb2gray(Im1);
Im2 = imread('thermalImage4.jpg');
Im2 = rgb2gray(Im2);

% Filtering
Im1_Gauss = imgaussfilt(Im1); % Gaussian
Im2_Gauss = imgaussfilt(Im2); % Gaussian

figure(1)
subplot(2, 1, 1); montage({Im1, Im1_Gauss}, 'size', [1 2]); xlabel('Original, Gaussian');
subplot(2, 1, 2); montage({Im2, Im2_Gauss}, 'size', [1 2]);

% =========================================================================
% Boundary Extraction & Region Filling
% =========================================================================
% Image 1
BW = imbinarize(Im1_Gauss);
[B, L] = bwboundaries(BW, 'noholes');

figure(2)
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

% Image 2
BW2 = imbinarize(Im2_Gauss);
[B2, L2] = bwboundaries(BW2, 'noholes');

figure(3)
imshow(label2rgb(L2, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B2)
    boundary2 = B2{k};
    plot(boundary2(:,2), boundary2(:,1), 'w', 'LineWidth', 2)
end

% =========================================================================
% Enhance Images
% =========================================================================
% Enhanced Image = Original + TopHat - BottomHat
seT = strel('disk', 5);
seT2 = strel('rectangle', [9 9]);
seT3 = strel('square', 7);
seT4 = strel('cube', 7);

tophat = imtophat(Im1_Gauss,seT);
tophat2 = imtophat(Im1_Gauss,seT2);
tophat3 = imtophat(Im1_Gauss,seT3);
tophat4 = imtophat(Im1_Gauss,seT4);

bottHat = imbothat(Im1_Gauss,seT);
bottHat2 = imbothat(Im1_Gauss,seT2);
bottHat3 = imbothat(Im1_Gauss,seT3);
bottHat4 = imbothat(Im1_Gauss,seT4);

enhanced = imsubtract(imadd(Im1_Gauss,tophat),bottHat);
enhanced2 = imsubtract(imadd(Im1_Gauss,tophat2),bottHat2);
enhanced3 = imsubtract(imadd(Im1_Gauss,tophat3),bottHat3);
enhanced4 = imsubtract(imadd(Im1_Gauss,tophat4),bottHat4);

figure(12)
subplot(4, 1, 1); montage({Im1_Gauss, tophat, bottHat, enhanced}, 'size', [1 4]);
subplot(4, 1, 2); montage({Im1_Gauss, tophat2, bottHat2, enhanced2}, 'size', [1 4]);
subplot(4, 1, 3); montage({Im1_Gauss, tophat3, bottHat3, enhanced3}, 'size', [1 4]);
subplot(4, 1, 4); montage({Im1_Gauss, tophat4, bottHat4, enhanced4}, 'size', [1 4]);

%{
% --------------
% Step 1: Get an input binary image of size NxN
% Im1_Gauss
% Step 2: Add impulse noise to the input binary image
Im1_Gauss_noise = imnoise(Im1_Gauss, 'salt & pepper', 0.4);
figure(13)
subplot(1, 1, 1); imshow(Im1_Gauss_noise); xlabel('Noisy Image');

% Step 3: Apply dilation to the noisy image
Im1Gauss_Noise = im2double(Im1_Gauss_noise);
se = strel('disk', 5.0);
Im1GaussNoiseD = imdilate(Im1Gauss_Noise, se);
edge = Im1GaussNoiseD - Im1Gauss_Noise;
BW = edge > 0.30;
figure(14)
subplot(1, 1, 1); montage({Im1Gauss_Noise, Im1GaussNoiseD, edge, BW}, 'size', [1 4]); 
xlabel('Original, Dilated, Edge=Dilated-Original, Edge>Threshold');

% Step 4: Apply erosion the resultant image
A0 = BW;
ser = strel('disk', 7.0); % Structuring element
Im1_erosion = imerode(A0, ser); % Erode the image by structuring element
figure(15)
subplot(2 ,1, 1); imshow(Im1_erosion); xlabel('Eroded Image');
subplot(2, 1, 2); imshow(A0-Im1_erosion); xlabel('Difference between binary image and erode image');

% Step 5: Apply erosion first time after closing
se0 = strel('disk', 7);
closeBW0 = imclose(Im1_erosion, se0);
figure(16)
imshow(closeBW0);
imshow(Im1_Gauss-closeBW6); % Gradient

% Step 6: Apply erosion two times after closing


% Step 7: Apply erosion three times after closing and obtain the noise free
% image


% Step 8: Repear step 2-7 by adding different impulse noise intensities
%}

% =========================================================================
% Detect Edges (Morphological Gradient)
% =========================================================================

%  ----
%  A. Using Greyscale and binary morphological operations, such as erosion,
%  dilation, opening and closing
%  ----

% EROSION
level = graythresh(Im1_Gauss);
% erosion with cross structuring element
length = 15;
NHOOD = zeros(length);
NHOOD(ceil(length/2), :) = 1;
NHOOD(:, ceil(length/2)) = 1;
se = strel('arbitrary', NHOOD);
erodedImg1 = imerode(Im1_Gauss, se);
% thresholding
BW = erodedImg1 > level;

figure(4)
subplot(3, 1, 1); imshow(Im1_Gauss); xlabel('Original');
subplot(3, 1, 2); imshow(erodedImg1); xlabel('Eroded Image');
subplot(3, 1, 3); imshow(uint8(BW)); xlabel('Eroded Image - Greythresh Level');

% ---------------------
A = Im1_Gauss;
s = strel('disk', 3.0); % Structuring element
C = imerode(A, s); % Erode the image by structuring element

A2 = Im1_Gauss;
s2 = strel('rectangle', [6 6]); % Structuring element
C2 = imerode(A2, s2); % Erode the image by structuring element

A3 = Im1_Gauss;
s3 = strel('square', 5.0); % Structuring element
C3 = imerode(A3, s3); % Erode the image by structuring element

A4 = Im1_Gauss;
s4 = strel('cube', 5.0); % Structuring element
C4 = imerode(A4, s4); % Erode the image by structuring element

A5 = Im1_Gauss;
s5 = strel('sphere', 4.0); % Structuring element
C5 = imerode(A5, s5); % Erode the image by structuring element

figure(5)
subplot(3, 3, 1); imshowpair(A, C); xlabel('Disk');
subplot(3, 3, 2); imshowpair(A2, C2); xlabel('Rectangle'); title('Highlighted Edge')
subplot(3, 3, 3); imshowpair(A3, C3); xlabel('Square');
subplot(3, 3, 4); imshowpair(A4, C4); xlabel('Cube');
subplot(3, 3, 5); imshowpair(A5, C5); xlabel('Sphere');

figure(6)
subplot(3, 3, 1); imshow(C); xlabel('Disk');
subplot(3, 3, 2); imshow(C2); xlabel('Rectangle'); title('Eroded Image');
subplot(3, 3, 3); imshow(C3); xlabel('Square');
subplot(3, 3, 4); imshow(C4); xlabel('Cube');
subplot(3, 3, 5); imshow(C5); xlabel('Sphere');

figure(30)
montage({C, C2, C3, C4, C5}, 'size', [1 5]); xlabel('Eroded Image');


% Difference between binary image and erode image
figure(7)
montage({(A-C), (A2-C2), (A3-C3), (A4-C4), (A5-C5)}, 'size', [1 5]); xlabel('Difference between binary image and erode image');

% ----------------------

% DILATION
% ref: Power point & scipt by Qiyuan Tian and David Chen
% change strel
Im1_Gauss = im2double(Im1_Gauss);
% dilation
se = strel('disk', 7.0);
Im1GaussD = imdilate(Im1_Gauss, se);

% subtraction and thresholding
edge = Im1GaussD - Im1_Gauss;
BW = edge > 0.15;

figure(8)
subplot(1, 1, 1); montage({Im1GaussD, edge, BW}, 'size', [1 3]); 
xlabel('Original, Dilated, Edge=Dilated-Original, Edge>Threshold');

% OPENING
% s = strel('rectangle', [40 30]);
s6 = strel('disk', 3);
s7 = strel('rectangle', [9 9]);
s8 = strel('square', 7.0);
s9 = strel('cube', 7.0);

BW6 = imopen(Im1_Gauss, s6);
BW7 = imopen(Im1_Gauss, s7);
BW8 = imopen(Im1_Gauss, s8);
BW9 = imopen(Im1_Gauss, s9);

figure(9)
%montage({BW6, BW7, BW8, BW9}, 'size', [1 4]);
% Gradient
montage({(Im1_Gauss-BW6), (Im1_Gauss-BW7), (Im1_Gauss-BW8), (Im1_Gauss-BW9)}, 'size', [1 4]);

% CLOSING
se6 = strel('disk', 3);
se7 = strel('rectangle', [9 9]);
se8 = strel('square', 7.0);
se9 = strel('cube', 7.0);

closeBW6 = imclose(Im1_Gauss, se6);
closeBW7 = imclose(Im1_Gauss, se7);
closeBW8 = imclose(Im1_Gauss, se8);
closeBW9 = imclose(Im1_Gauss, se9);

figure(10)
%montage({closeBW6, closeBW7, closeBW8, closeBW9}, 'size', [1 4]);
% Gradient
montage({(Im1_Gauss-closeBW6), (Im1_Gauss-closeBW7), (Im1_Gauss-closeBW8), (Im1_Gauss-closeBW9)}, 'size', [1 4]);
%  ----
%  B. Different structural elements (Multiscale Edge Detectors)
%  ----
% shown on power point

%  ----
%  C. A comparison of above morphological techniques in the taxonomy of your
%  applications
%  ----
% on power point