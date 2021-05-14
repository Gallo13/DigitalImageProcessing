% Created by: Jessica Gallo
% Date created: 5/7/2021
% Last modified: 5/13/2021
% CSC741 Digital Image Processing
% Final Project
% Face Detection of different emotions with thermal images

% -----------
% To Do:
%   - load all folders into dataset 
%   - filter all images 
%   - edge detection and morphological operations 
%   - save images to .TIFF
%   - save new filtered dataset to local and use in CNN
% ------------

clear;
clc;
close all;

% =========================================================================
% Load Dataset, Gradient Filter, Histogram Equalization, Sobel Edge
% Detector, JPG to TIFF
% =========================================================================
%{
imds = imageDatastore('C:\Users\User\Desktop\TD_IR_RGB_CROPPED\RGB-faces-128x128\', ...
    'IncludeSubFolders', true, 'FileExtensions', '.jpg', 'LabelSource', 'foldernames');
figure(1)
imshow(preview(imds));

length = length(imds);
%length = imds.Count;

for i = 1:length
    img = readimage(imds, i);
    filtered = medfilt2(img, [3 3]);
end
%}

% Load 'Neutral' Images
% ---------------------
D = 'C://Users//User//Desktop//Thermal//1_Neutral';
feature_matrix_Neutral = [];
Files=dir('C://Users//User//Desktop//Thermal//1_Neutral');
for k=3:length(Files)
   
   Files(k);
   FileNames=Files(k).name;
   S = fullfile(D,FileNames); 
   I = imread(S);
   
   %figure(1)
   %imshow(I);
   %drawnow;
   
   % Gaussian Filter to Images (maybe change to Median Filter)
   % ---------------------------------------------------------
   f = imgaussfilt(I);
   
   %figure(2)
   %imshow(f);
   %drawnow;
   
   % Histogram Equalization
   % ----------------------
   fhist = histeq(f);
   
   %figure(3)
   %imshow(fhist);
   %drawnow;
   
   % Sobel Edge Detection
   % --------------------
   % RGB to Grey
   fhistG = fhist(:,:,1);
   
   %figure(4)
   %imshow(fhistG);
   %drawnow;
   
   A = fhistG;

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

   sobel = uint8(Im1_Sobel);
   
   %figure(5)
   %imshow(sobel);
   %drawnow;
   
   %feature_matrix_Neutral = [feature_matrix_Neutral;sobel];
   
   %figure(6)
   %imshow(feature_matrix_Neutral);
   %drawnow;
   
   % Save the processed images 
   %save(FileNames,'sobel');
   imwrite(sobel, FileNames);
end
% Check Images
%{
figure(1)
subplot(7, 1, 1); imshow(S); xlabel('S');
subplot(7, 1, 2); imshow(I); xlabel('I');
subplot(7, 1, 7); imshow(I); drawnow; xlabel('I');
subplot(7, 1, 3); imshow(f); xlabel('Filtered Im');
subplot(7, 1, 4); imshow(fhist); xlabel('HistEq Im');
subplot(7, 1, 5); imshow(fhistG); xlabel('HistEq Im Grey');
subplot(7, 1, 6); imshow(sobel); xlabel('Sobel');
%}

% Load 'Smiling' Images
% ---------------------
D = 'C://Users//User//Desktop//Thermal//2_Smiling';
feature_matrix_Smiling = [];
Files=dir('C://Users//User//Desktop//Thermal//2_Smiling');
for k=3:length(Files)
   
   FileNames=Files(k).name;
   S = fullfile(D,FileNames);
   I = imread(S);
   
   % Gaussian Filter to Images (maybe change to Median Filter)
   % ---------------------------------------------------------
   f = imgaussfilt(I);
   
   % Histogram Equalization
   % ----------------------
   fhist = histeq(f);
   
   % Sobel Edge Detection
   % --------------------
   % RGB to Grey
   fhistG = fhist(:,:,1);
   A = fhistG;

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

   sobel = uint8(Im1_Sobel);
   
   %feature_matrix_Smiling = [feature_matrix_Smiling;sobel];
   
   % Save the processed images 
   imwrite(sobel, FileNames);
end
% Check Images
%{
figure(21)
subplot(6, 1, 1); imshow(S); xlabel('S');
subplot(6, 1, 2); imshow(I); xlabel('I');
subplot(6, 1, 3); imshow(f); xlabel('Filtered Im');
subplot(6, 1, 4); imshow(fhist); xlabel('HistEq Im');
subplot(6, 1, 5); imshow(fhistG); xlabel('HistEq Im Grey');
subplot(6, 1, 6); imshow(sobel); xlabel('Sobel');
%}
% Load 'Sleepy' Images
% --------------------
D = 'C://Users//User//Desktop//Thermal//3_Sleepy';
feature_matrix_Sleepy = [];
Files=dir('C://Users//User/Desktop//Thermal//3_Sleepy');
for k=3:length(Files)
   
   FileNames=Files(k).name;
   S = fullfile(D,FileNames);
   I = imread(S);
   
   % Gaussian Filter to Images (maybe change to Median Filter)
   % ---------------------------------------------------------
   f = imgaussfilt(I);
   
   % Histogram Equalization
   % ----------------------
   fhist = histeq(f);
   
   % Sobel Edge Detection
   % --------------------
   % RGB to Grey
   fhistG = fhist(:,:,1);
   A = fhistG;

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

   sobel = uint8(Im1_Sobel);
   
   %feature_matrix_Sleepy = [feature_matrix_Sleepy;sobel];
   
   % Save the processed images 
   imwrite(sobel, FileNames);
   
end
%{
% Check Images
figure(8)
%subplot(6, 1, 1); imshow(S); xlabel('S');
%subplot(6, 1, 2); imshow(I); xlabel('I');
%subplot(6, 1, 3); imshow(f); xlabel('Filtered Im');
%subplot(6, 1, 4); imshow(fhist); xlabel('HistEq Im');
subplot(1, 1, 1); imshow(fhistG); xlabel('HistEq Im Grey');
%subplot(6, 1, 6); imshow(sobel); xlabel('Sobel');
%}
% Load 'Surprise' Images
% ----------------------
D = 'C://Users//User//Desktop//Thermal//4_Surprise';
feature_matrix_Surprise = [];
Files=dir('C://Users//User//Desktop//Thermal//4_Surprise');
for k=3:length(Files)
   
   FileNames=Files(k).name;
   S = fullfile(D,FileNames);
   I = imread(S);
   
   % Gaussian Filter to Images (maybe change to Median Filter)
   % ---------------------------------------------------------
   f = imgaussfilt(I);
   
   % Histogram Equalization
   % ----------------------
   fhist = histeq(f);
   
   % Sobel Edge Detection
   % -------------------- 
   % RGB to Grey
   fhistG = fhist(:,:,1);
   A = fhistG;

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

   sobel = uint8(Im1_Sobel);
   
   %feature_matrix_Surprise = [feature_matrix_Surprise;sobel];
   
   % Save the processed images 
   imwrite(sobel, FileNames);
   
end
%{
% Check Images
figure(14)
subplot(6, 1, 1); imshow(S); xlabel('S');
subplot(6, 1, 2); imshow(I); xlabel('I');
subplot(6, 1, 3); imshow(f); xlabel('Filtered Im');
subplot(6, 1, 4); imshow(fhist); xlabel('HistEq Im');
subplot(6, 1, 5); imshow(fhistG); xlabel('HistEq Im Grey');
subplot(6, 1, 6); imshow(sobel); xlabel('Sobel');
%}
% Load 'Sunglasses' Images
% ------------------------
D = 'C://Users//User//Desktop//Thermal//5_Sunglasses';
feature_matrix_Sunglasses = [];
Files=dir('C://Users//User//Desktop//Thermal//5_Sunglasses');
for k=3:length(Files)
   
   FileNames=Files(k).name;
   S = fullfile(D,FileNames);
   I = imread(S);
   
   % Gaussian Filter to Images (maybe change to Median Filter)
   % ---------------------------------------------------------
   f = imgaussfilt(I);
   
   % Histogram Equalization
   % ----------------------
   fhist = histeq(f);
   
   % Sobel Edge Detection
   % --------------------
   % RGB to Grey
   fhistG = fhist(:,:,1);
   A = fhistG;

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

   sobel = uint8(Im1_Sobel);
   
   %feature_matrix_Sunglasses = [feature_matrix_Sunglasses;sobel];
   
   % Save the processed images 
   imwrite(sobel, FileNames);
end
%{
% Check Images
figure(51)
subplot(6, 1, 1); imshow(S); xlabel('S');
subplot(6, 1, 2); imshow(I); xlabel('I');
subplot(6, 1, 3); imshow(f); xlabel('Filtered Im');
subplot(6, 1, 4); imshow(fhist); xlabel('HistEq Im');
subplot(6, 1, 5); imshow(fhistG); xlabel('HistEq Im Grey');
subplot(6, 1, 6); imshow(sobel); xlabel('Sobel');
%}