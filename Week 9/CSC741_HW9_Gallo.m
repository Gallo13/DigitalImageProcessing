% Created by: Jessica Gallo
% Date Created: 4/12/2021
% Last Modified: 4/14/2021
% Homework 9
% Image Segmentation Assignment

clc;
close all;

% =========================================================================
% Implement histogram-based segmentation on your image. Identify the peaks
% of your histogram with the "objects" that they correspond to. Show your
% image, its histogram, the ranges, etc. Show the identified objects.
% =========================================================================
Im = imread("thermalImage3.jpg");
Im = rgb2gray(Im);

Im2 = imread("00483_s_20aqapbvgk0483.jpg");
Im2 = rgb2gray(Im2);

figure(1)
subplot(2, 2, 1); imshow(Im2);
subplot(2, 2, 2); imhist(Im2);
subplot(2, 2, 3); imshow(Im);
subplot(2, 2, 4); imhist(Im);

% =========================================================================
% Implement a global thresholding algorithm and experiment on images with
% (a) two major modes (b) three significant modes and observe the results
% using different initial estimates of T. Discuss your findings
% =========================================================================
% Histogram Based Segmentation with Otsu's Method |
% -------------------------------------------------

% Calculate a 255-bin histogram for the image
[counts, x] = imhist(Im, 255);
figure(2)
stem(x, counts);

% 2 major modes with otsuthresh
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute a global threshold using the histogram counts
T1 = otsuthresh(counts);
% Create a binary image using the computed threshold
BW1 = imbinarize(Im, T1);

% Display the image
figure(3); 
subplot(1, 1, 1); imshowpair(Im, BW1, 'montage'); xlabel('Original Image & Segmented Image with Otsu Thresholding');

% 2 major modes with greythresh
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute a global threshold using the histogram counts
T2 = graythresh(counts); % T2 = 0.2
% Create a binary image using the computed threshold
BW2 = imbinarize(Im, T2);

% Display the image
figure(4); 
subplot(1, 1, 1); imshowpair(Im, BW2, 'montage'); xlabel('Original Image & Segmented Image with Grey Thresholding');

% 2 major modes with imbinarize
% Iterative Global Thresholding
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
T3=0.5*(double(min(Im(:)))+ double(max(Im(:))));
done = false;
while ~done
g= Im>=T3;
Tnext=0.5*(mean(Im(g)) + mean(Im(~g)));
done=abs(T3-Tnext)>0.5;
T3=Tnext;
end

BW3 = imbinarize(Im, T3);
figure(5);
subplot(1, 1, 1); imshowpair(Im, BW3, 'montage'); xlabel('Original & Segmented Images with Iterative Global Thresholding');

% Summary
%figure(6)
%subplot(1, 1, 1); montage(Im, BW1, BW2, BW3, 'size', [1 4]); xlabel ('Original, Otsu, Grey, Iterative Global Thresholding');

% -----------------------------------------------------------
% 3 major modes with greythresh
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Calculate a 255-bin histogram for the image
[counts, x] = imhist(Im2, 255);
%figure(7)
%stem(x, counts);

% 3 major modes with otsuthresh
% ~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute a global threshold using the histogram counts
T4 = otsuthresh(counts);
% Create a binary image using the computed threshold
BW4 = imbinarize(Im2, T4);

% Display the image
figure(8); 
subplot(1, 1, 1); imshowpair(Im2, BW4, 'montage'); xlabel('Original Image & Segmented Image with Otsu Thresholding');

% 3 major modes with greythresh
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute a global threshold using the histogram counts
T5 = graythresh(counts); % T2 = 0.2
% Create a binary image using the computed threshold
BW5 = imbinarize(Im2, T5);

% Display the image
figure(9); 
subplot(1, 1, 1); imshowpair(Im2, BW5, 'montage'); xlabel('Original Image & Segmented Image with Grey Thresholding');

% 3 major modes with imbinarize
% Iterative Global Thresholding
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
T6=0.5*(double(min(Im2(:)))+ double(max(Im2(:))));
done = false;
while ~done
g= Im2>=T6;
Tnext=0.5*(mean(Im2(g)) + mean(Im2(~g)));
done=abs(T6-Tnext)>0.5;
T6=Tnext;
end

BW6 = imbinarize(Im2, T1);
figure(10);
subplot(1, 1, 1); imshowpair(Im2, BW6, 'montage'); xlabel('Original & Segmented Images with Iterative Global Thresholding');

% Summary
figure(11)
subplot(1, 1, 1); montage(Im2, BW4, BW5, BW6, 'size', [1 4]); xlabel ('Original, Otsu, Grey, Iterative Global Thresholding');


% =========================================================================
% Implement multilevel and adaptive thresholding algorithms and experiment
% on images and observe the results using different T's initial estimates.
% =========================================================================
% Multilevel Thresholding (2 Regions)
level = multithresh(Im, 4);
seg_Im = imquantize(Im,level);
figure(12)
subplot(1, 1, 1); imshow(seg_Im,[]); xlabel('Multilevel Thresholding into 2 Regions');

level2 = multithresh(Im2, 4);
seg_Im2 = imquantize(Im2,level2);
figure(13)
subplot(1, 1, 1); imshow(seg_Im2,[]); xlabel('Multilevel Thresholding into 2 Regions');

% Multilevel Thresholding (3 Levels/2 Thresholds)
thresh = multithresh(Im, 6);
seg_Im3 = imquantize(Im,thresh);
RGB = label2rgb(seg_Im3);

thresh2 = multithresh(Im2, 6);
seg_Im4 = imquantize(Im2,thresh2);
RGB2 = label2rgb(seg_Im4);

%figure(14)
%subplot(1, 1, 1); imshow(RGB2); xlabel('Multilevel Thresholding into 2 Levels');

% Adaptive Thresholding
T7 = adaptthresh(Im, 0.1);
BW7 = imbinarize(Im,T7);

T8 = adaptthresh(Im2, 0.1);
BW8 = imbinarize(Im2,T8);

T9 = adaptthresh(Im, 0.4);
BW9 = imbinarize(Im,T9);

T10 = adaptthresh(Im2, 0.4);
BW10 = imbinarize(Im2,T10);

T11 = adaptthresh(Im, 0.6);
BW11 = imbinarize(Im,T11);

T12 = adaptthresh(Im2, 0.6);
BW12 = imbinarize(Im2,T12);

figure(14)
subplot(2, 1, 1); montage({BW7, BW9, BW11}, 'Size', [1 3]); xlabel('Adaptive Thresholding');
figure(15)
subplot(2, 1, 1); montage({BW8, BW10, BW12}, 'Size', [1 3]); xlabel('Adaptive Thresholding');

% BONUSES ----------------------------------------------------------------
% =========================================================================
% Implement mean histogram stretching for image enhancement and sementation
% =========================================================================
% noise reduction
ImBi = imbinarize(Im);
% Median Filter
Imfilter = medfilt2(Im);
% contrast reduction
Imhisteq = histeq(Imfilter);

% Adaptive Binarization
ImhisteqBi = imbinarize(Imhisteq, 'adaptive', 'Sensitivity', 0.5);
figure(16)
subplot(1, 1, 1); montage({Im, ImBi, Imfilter, Imhisteq, ImhisteqBi}, 'Size', [1 5]); xlabel('Original, Binary, Medium Filter, Histogram Stretch, Binary Enhanced (Filter/HistStretch)')

ImhisteqBi3 = imbinarize(Imhisteq, 'adaptive', 'Sensitivity', 0.1);
ImhisteqBi4 = imbinarize(Imhisteq, 'adaptive', 'Sensitivity', 0.4);
ImhisteqBi5 = imbinarize(Imhisteq, 'adaptive', 'Sensitivity', 0.6);

figure(17)
subplot(1, 1, 1); montage({Im, ImhisteqBi3, ImhisteqBi4, ImhisteqBi5}, 'Size', [1 4]); xlabel('Original, 0.1, 0.4, 0.6');

% noise reduction
ImBi2 = imbinarize(Im2);
% Median Filter
Imfilter2 = medfilt2(Im2);
% contrast reduction
Imhisteq2 = histeq(Imfilter2);
% Adaptive Binarization
ImhisteqBi2 = imbinarize(Imhisteq2, 'adaptive', 'Sensitivity', 0.5);
figure(18)
subplot(1, 1, 1); montage({Im2, ImBi2, Imfilter2, Imhisteq2, ImhisteqBi2}, 'Size', [1 5]); xlabel('Original, Binary, Medium Filter, Histogram Stretch, Binary Enhanced (Filter/HistStretch)')

ImhisteqBi6 = imbinarize(Imhisteq2, 'adaptive', 'Sensitivity', 0.1);
ImhisteqBi7 = imbinarize(Imhisteq2, 'adaptive', 'Sensitivity', 0.4);
ImhisteqBi8 = imbinarize(Imhisteq2, 'adaptive', 'Sensitivity', 0.6);

figure(19)
subplot(1, 1, 1); montage({Im2, ImhisteqBi6, ImhisteqBi7, ImhisteqBi8}, 'Size', [1 4]); xlabel('Original, 0.1, 0.4, 0.6');

% =========================================================================
% Implement adaptive gaussian thresholding algorithm
% =========================================================================
% MAIN TO EXPERIMENT WITH
% 0.4/bright
T13 = adaptthresh(Im, 0.1, 'ForegroundPolarity', 'bright', 'NeighborhoodSize', 2*floor(size(Im)/16)+1, "Statistic", 'gaussian');
BW13 = imbinarize(Im, T13);
% 0.4 dark
T14 = adaptthresh(Im2, 0.1, 'ForegroundPolarity','bright', 'NeighborhoodSize', 2*floor(size(Im)/16)+1, "Statistic", 'gaussian');
BW14 = imbinarize(Im2, T14);

% 0.2/dark
T15 = adaptthresh(Im, 0.4, 'ForegroundPolarity','bright', 'NeighborhoodSize', 2*floor(size(Im)/16)+1, "Statistic", 'gaussian');
BW15 = imbinarize(Im, T15);
% 0.2/bright
T16 = adaptthresh(Im2, 0.4, 'ForegroundPolarity','bright', 'NeighborhoodSize', 2*floor(size(Im)/16)+1, "Statistic", 'gaussian');
BW16 = imbinarize(Im2, T16);

% 0.7/dark
T17 = adaptthresh(Im, 0.6, 'ForegroundPolarity','bright', 'NeighborhoodSize', 2*floor(size(Im)/16)+1, "Statistic", 'gaussian');
BW17 = imbinarize(Im, T17);
% 0.7/bright
T18 = adaptthresh(Im2, 0.6, 'ForegroundPolarity','bright', 'NeighborhoodSize', 2*floor(size(Im)/16)+1, "Statistic", 'gaussian');
BW18 = imbinarize(Im2, T18);

figure(31)
subplot(1, 1, 1); montage({BW13, BW15, BW17}, 'Size', [1 3]);

figure(32)
subplot(1, 1, 1); montage({BW14, BW16, BW18}, 'Size', [1 3]);


