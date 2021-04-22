% Created by: Jessica Gallo
% Date created: 4/16/2021
% Last modified: 4/22/2021
% CSC741 Digital Image Processing
% Project 1
% Edge Detection

clc;
close all;

Im1 = imread('thermalImage3.jpg');
Im1 = rgb2gray(Im1);
Im2 = imread('00483_s_20aqapbvgk0483.jpg');
Im2 = rgb2gray(Im2);
Im2 = imresize(Im2, [128, 128]);

figure(1)
subplot(2, 1, 1); imshow(Im1);
subplot(2, 1, 2); imshow(Im2);

% Filtering and Enhancing ------------------------------------------------

% Gaussian Filter
Im1_denoise_G = imgaussfilt(Im1, 0.5);
Im2_denoise_G = imgaussfilt(Im2, 0.5);

% Median Mean Filter
% original
Im1_denoise_Med = medfilt2(Im1);
Im2_denoise_Med = medfilt2(Im2);

% Contrast Enhancement
% Histogram Equalization/Stretching
ImHistEq1 = histeq(Im1);
ImHistEqMed1 = histeq(Im1_denoise_Med);

ImHistEq2 = histeq(Im2);
ImHistEqMed2 = histeq(Im2_denoise_Med);

ImHistEqG1 = histeq(Im1_denoise_G);
ImHistEqG2 = histeq(Im2_denoise_G);

figure(2)
subplot(4, 1, 1); montage({Im1, Im1_denoise_G, Im1_denoise_Med}, 'size', [1 3]); xlabel('Original, Gaussian, Median');
subplot(4, 1, 2); montage({ImHistEq1, ImHistEqMed1, ImHistEqG1}, 'size', [1 3]); xlabel('Original HistEq, Gaussian HistEq, Median HistEq');
subplot(4, 1, 3); montage({Im2, Im2_denoise_G, Im2_denoise_Med}, 'size', [1 3]); xlabel('Original, Gaussian, Median');
subplot(4, 1, 4); montage({ImHistEq2, ImHistEqMed2, ImHistEqG2}, 'size', [1 3]); xlabel('Original HistEq, Gaussian HistEq, Median HistEq');

% A -----------------------------------------------------------------------
% =========================================================================
% Sobel edge detector with threshold = 20,40,128. Apply these to various
% visible and thermal images
% =========================================================================
% ref: https://www.imageeprocessing.com/2013/07/sobel-edge-detection-part-2.html
% without using edge fuction

% Thermal Images -------------------------
% re-run code and change variable A to images you want & change threshold
% amount. Everything else stays the same except vairables A and threshold
A = ImHistEqG1;

Im1_Sobel = zeros(size(A));

F1 = [-1 0 1; -2 0 2; -1 0 1];
F2 = [-1 -2 -1; 0 0 0; 1 2 1];

A = double(A);

for i = 1:size(A, 1)-2
    for j = 1:size(A, 2)-2
        % Gradient Operations
        Gx = sum(sum(F1.*A(i:i+2, j:j+2)));
        Gy = sum(sum(F2.*A(i:i+2, j:j+2)));
        
        % Magnitude of Vector
        Im1_Sobel(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
    end
end

I = uint8(Im1_Sobel);

% Define Threhold value
% Thermal / T=20
threshold = 20;
B = max(Im1_Sobel, threshold);
B(B==round(threshold))=0;

B = imbinarize(B);
figure(3), imshow(B); title('Sobel Gaussian Thermal Image HistEq/T=20');

% Visual Images -------------------------
% re-run code and change variable A to images you want & change threshold
% amount. Everything else stays the same except vairables A and threshold
A2 = ImHistEqG2;

Im2_Sobel = zeros(size(A2));

F1 = [-1 0 1; -2 0 2; -1 0 1];
F2 = [-1 -2 -1; 0 0 0; 1 2 1];

A2 = double(A2);

for i = 1:size(A2, 1)-2
    for j = 1:size(A2, 2)-2
        % Gradient Operations
        Gx = sum(sum(F1.*A2(i:i+2, j:j+2)));
        Gy = sum(sum(F2.*A2(i:i+2, j:j+2)));
        
        % Magnitude of Vector
        Im2_Sobel(i+1, j+1) = sqrt(Gx.^2 + Gy.^2);
    end
end

I2 = uint8(Im2_Sobel);

% Define Threhold value
% Thermal / T=20
threshold2 = 20;
B2 = max(Im2_Sobel, threshold2);
B2(B2==round(threshold2))=0;

B2 = imbinarize(B2);
figure(4), imshow(B2); title('Sobel Visual Gaussian Image HistEq/T=20');


% =========================================================================
% Zero-crossing on the following four types of images to get edge images
% (choose proper thresholds)
% =========================================================================

powerline1 = imread('powerline1.png');
powerline1 = rgb2gray(powerline1);
powerline1G = imgaussfilt(powerline1); % gaussian filter

powerline2 = imread('powerline2.png');
powerline2 = rgb2gray(powerline2);
powerline2G = imgaussfilt(powerline2);

powerline3 = imread('powerline3.png');
powerline3 = rgb2gray(powerline3);
powerline3G = imgaussfilt(powerline3);

powerline4 = imread('powerline4.png');
powerline4 = rgb2gray(powerline4);
powerline4G = imgaussfilt(powerline4);

% figure out what threshold is

BW1 = edge(powerline1, 'zerocross', 0);
BW2 = edge(powerline2, 'zerocross', 0);
BW3 = edge(powerline3, 'zerocross', 0);
BW4 = edge(powerline4, 'zerocross', 0);
BW5 = edge(powerline1, 'zerocross', 0.01);
BW6 = edge(powerline2, 'zerocross', 0.01);
BW7 = edge(powerline3, 'zerocross', 0.01);
BW8 = edge(powerline4, 'zerocross', 0.01);

figure(5)
subplot(4, 1, 1); montage({powerline1, powerline1G, BW1, BW5}, 'size', [1 4]);
subplot(4, 1, 2); montage({powerline2, powerline2G, BW2, BW6}, 'size', [1 4]);
subplot(4, 1, 3); montage({powerline3, powerline3G, BW3, BW7}, 'size', [1 4]);
subplot(4, 1, 4); montage({powerline4, powerline4G, BW4, BW8}, 'size', [1 4]);

% =========================================================================
% Canny edge detector
% ========================================================================
BW9 = edge(Im1, 'canny', 0);
BW10 = edge(Im1_denoise_G, 'canny', 0.4);
BW11 = edge(ImHistEq1, 'canny', 0.4);
BW12 = edge(ImHistEqG1, 'canny', 0.4);

BW13 = edge(Im2, 'canny', 0.4);
BW14 = edge(Im2_denoise_G, 'canny', 0.4);
BW15 = edge(ImHistEq2, 'canny', 0.4);
BW16 = edge(ImHistEqG2, 'canny', 0.4);

figure(6)
subplot(4, 1, 1); montage({BW9, BW10, BW11, BW12}, 'size', [1 4]);
subplot(4, 1, 3); montage({BW13, BW14, BW15, BW16}, 'size', [1 4]);

% =========================================================================
% Robinson edge detector
% ========================================================================
%{
direction_list = ["North", "West", "South", "East", ...
    "NorthWest", "SouthWest", "SouthEast", "NorthEast"];

IM=double(Im1);
g1=[-1 0 1; -2 0 2; -1 0 1]; %N
g2=[0 1 2; -1 0 1; -2 -1 0]; %W
g3=[1 2 1; 0 0 0; -1 -2 -1]; %S
g4=[2 1 0; 1 0 -1; 0 -1 -2]; %E
g5=[1 0 -1; 2 0 -2; 1 0 -1]; %NW
g6=[0 -1 -2; 1 0 -1;2 1 0]; %SW
g7=[-1 -2 -1; 0 0 0;1 2 1]; %SE
g8=[-2,-1,0; -1,0,1;0,1,2]; % NE
IM1_h(:,:,1)=imfilter(IM,g1,'replicate');
IM1_h(:,:,2)=imfilter(IM,g2,'replicate');
IM1_h(:,:,3)=imfilter(IM,g3,'replicate');
IM1_h(:,:,4)=imfilter(IM,g4,'replicate');
IM1_h(:,:,5)=imfilter(IM,g5,'replicate');
IM1_h(:,:,6)=imfilter(IM,g6,'replicate');
IM1_h(:,:,7)=imfilter(IM,g7,'replicate');
IM1_h(:,:,8)=imfilter(IM,g8,'replicate');


figure(7)
for i=1:8
    subplot(1, 1, 1); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end

y1=max(IM1_h(:,:,1),IM1_h(:,:,2));
y2=max(y1,IM1_h(:,:,3));
y3=max(y2,IM1_h(:,:,4));
y4=max(y3,IM1_h(:,:,5));
y5=max(y4,IM1_h(:,:,6));
y6=max(y5,IM1_h(:,:,7));
y7=max(y6,IM1_h(:,:,8));
y=y7;

figure(8)
subplot(2, 1, 1); imshow(uint8(y));title('Robinson Operator Image')
j = 2;
for T = 30:30:90
    subplot(2, 1, 2);imshow(uint8(threshold_IM1(y,T)))
    title(['Binary Image, T = ' num2str(T)])
    j=j+1;
end
%}
% B -----------------------------------------------------------------------
% =========================================================================
% Frei and Chen gradient operator with threshold = 20
% =========================================================================
IM=double(Im1); % change this image
g1 = [1, sqrt(2), 1;0, 0, 0;-1, -sqrt(2), -1];
g2 = [1, 0, -1;sqrt(2), 0, -sqrt(2);1, 0, -1];
g3 = [0, -1, sqrt(2);1, 0, -1;-sqrt(2), 1, 0];
g4 = [sqrt(2), -1, 0;-1, 0, 1;0, 1, -sqrt(2)];
g5 = [0, 1, 0;-1, 0, -1;0, 1, 0];
g6 = [-1, 0, 1;0, 0, 0;1, 0, -1];
g7 = [1, -2, 1;-2, 4, -2;1, -2, 1];
g8 = [-2, 1, -2;1, 4, 1;-2, 1, -2];
g9 = ones(3,3)/9;
IM1_h(:,:,1)=imfilter(IM,g1,'replicate');
IM1_h(:,:,2)=imfilter(IM,g2,'replicate');
IM1_h(:,:,3)=imfilter(IM,g3,'replicate');
IM1_h(:,:,4)=imfilter(IM,g4,'replicate');
IM1_h(:,:,5)=imfilter(IM,g5,'replicate');
IM1_h(:,:,6)=imfilter(IM,g6,'replicate');
IM1_h(:,:,7)=imfilter(IM,g7,'replicate');
IM1_h(:,:,8)=imfilter(IM,g8,'replicate');
IM1_h(:,:,9)=imfilter(IM,g9,'replicate');

figure(7)
for i=1:9
subplot(3, 3, 1); imshow(uint8(IM1_h(:,:,i))); title(['M' num2str(i)]);
end

y1=max(IM1_h(:,:,1),IM1_h(:,:,2));
y2=max(y1,IM1_h(:,:,3));
y3=max(y2,IM1_h(:,:,4));
y4=max(y3,IM1_h(:,:,5));
y5=max(y4,IM1_h(:,:,6));
y6=max(y5,IM1_h(:,:,7));
y7=max(y6,IM1_h(:,:,8));
y8=max(y7,IM1_h(:,:,9));
y=y8;

figure(8)
subplot(2, 1, 1); imshow(uint8(y));title('Frei Chen Image')
j = 2;
for T = 100
subplot(2, 1, 2); imshow(uint8(threshold_IM2(y,T)))
title(['Binary Image, T = ' num2str(T)])
j=j+1;
end

% =========================================================================
% Kirsch compass operator with threshold = 20
% =========================================================================
direction_list = ["North", "West", "South", "East", ...
    "NorthWest", "SouthWest", "SouthEast", "NorthEast"];

IM=double(ImHistEqG2); % change image
g1=[5,5,5; -3,0,-3; -3,-3,-3];
g2=[5,5,-3; 5,0,-3; -3,-3,-3];
g3=[5,-3,-3; 5,0,-3; 5,-3,-3];
g4=[-3,-3,-3; 5,0,-3; 5,5,-3];
g5=[-3,-3,-3; -3,0,-3; 5,5,5];
g6=[-3,-3,-3; -3,0,5;-3,5,5];
g7=[-3,-3,5; -3,0,5;-3,-3,5];
g8=[-3,5,5; -3,0,5;-3,-3,-3];
IM1_h(:,:,1)=imfilter(IM,g1,'replicate');
IM1_h(:,:,2)=imfilter(IM,g2,'replicate');
IM1_h(:,:,3)=imfilter(IM,g3,'replicate');
IM1_h(:,:,4)=imfilter(IM,g4,'replicate');
IM1_h(:,:,5)=imfilter(IM,g5,'replicate');
IM1_h(:,:,6)=imfilter(IM,g6,'replicate');
IM1_h(:,:,7)=imfilter(IM,g7,'replicate');
IM1_h(:,:,8)=imfilter(IM,g8,'replicate');

figure(9)
for i=1
    subplot(4, 4, 1); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=2
    subplot(4, 4, 2); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=3
    subplot(4, 4, 3); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=4
   subplot(4, 4, 4); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=5
    subplot(4, 4, 5); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=6
    subplot(4, 4, 6); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=7
    subplot(4, 4, 7); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end
for i=8
   subplot(4, 4, 8); imshow(uint8(IM1_h(:,:,i))); title(direction_list(i));
end

%imshow(uint8(IM_kirsch))


y1=max(IM1_h(:,:,1),IM1_h(:,:,2));
y2=max(y1,IM1_h(:,:,3));
y3=max(y2,IM1_h(:,:,4));
y4=max(y3,IM1_h(:,:,5));
y5=max(y4,IM1_h(:,:,6));
y6=max(y5,IM1_h(:,:,7));
y7=max(y6,IM1_h(:,:,8));
y=y7;

%figure(2)
%imshow(uint8(y))
%title('Kirsch Image')


figure(10)
subplot(2, 1, 1); imshow(uint8(y));title('Kirsch Image')
j = 2;
for T = 20 % change threshold
subplot(2, 1, 2); imshow(uint8(threshold_IM3(y,T)))
title(['Binary Image, T = ' num2str(T)])
j=j+1;
end

% =========================================================================
% Functions
% =========================================================================
%{
% Robinson Compass Operator
function T_IM = threshold_IM1(Im,thresholdR)
[N, M] = size(Im1);
for i = 1:N
    for j = 1:M
    if (Im1(i,j) >= thresholdR)
        Im1(i,j) = 255;   
    else
        Im1(i,j) = 0;
    end
    end
end
T_IM = IM;
end
%}

% Frei Chen
function T_IM = threshold_IM2(IM, thresholdFC)
[N, M] = size(IM);
for i = 1:N
    for j = 1:M
    if (IM(i,j) >= thresholdFC)
        IM(i,j) = 255;   
    else
        IM(i,j) = 0;
    end
    end
end
T_IM = IM;
end

% Kirsch
function T_IM = threshold_IM3(IM,thresholdK)
[N, M] = size(IM);
for i = 1:N
    for j = 1:M
    if (IM(i,j) >= thresholdK)
        IM(i,j) = 255;   
    else
        IM(i,j) = 0;
    end
    end
end
T_IM = IM;
end