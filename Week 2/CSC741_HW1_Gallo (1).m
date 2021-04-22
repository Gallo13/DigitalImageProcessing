close all; clear;

% 2 images from Stanford Cars Image Dataset

% Defining Images
Im1 = imread('00001.jpg');
Im2 = imread('00002.jpg');

% Show color and grey image
figure(1)

subplot(2, 2, 1); imshow(Im1); title('Color Image 1');
Im1 = rgb2gray(Im1);
subplot(2, 2, 2); imshow(Im1); title('Greyscale Image 1');

subplot(2, 2, 3); imshow(Im2); title('Color Image 2');
Im2 = rgb2gray(Im2);
subplot(2, 2, 4); imshow(Im2); title('Greyscale Image 2');

% Image change Image Resolution 
% Low Res
num_col = 32;          % change resolution
num_row = 32;          % change resolution
I1ResizeLow = imresize(Im1,[num_col num_row]);
I2ResizeLow = imresize(Im2,[num_col num_row]);

figure(2)
subplot(2,1,1); imshow(I1ResizeLow); title([num2str(num_row) 'x' num2str(num_col)])
subplot(2,1,2); imshow(I2ResizeLow);

% High Res
num_col = 1028;          % change resolution
num_row = 1028;          % change resolution
I1ResizeHigh = imresize(Im1,[num_col num_row]); 
I2ResizeHigh = imresize(Im2,[num_col num_row]); 

figure(3)
subplot(2,1,1); imshow(I1ResizeHigh); title([num2str(num_row) 'x' num2str(num_col)]);
subplot(2,1,2); imshow(I2ResizeHigh);

% Image Scaling/Resizing
% images must be same size to add
Im1 = imresize(Im1, [500, 700]);
Im2 = imresize(Im2, [500, 700]);

% Image Addition
ImAdd = imadd(Im1,Im2);

weight1 = 0.2;
weight2 = 1 - weight1;

figure(4)
subplot(2, 2, 1); imshow(Im1); xlabel('Original Image 1');
title('Image Addition');
subplot(2, 2, 2); imshow(Im2); xlabel('Original Image 2');
subplot(2, 2, 3); imshow(ImAdd,[]); xlabel('Regular Addition');
subplot(2, 2, 4); imshow(weight1*Im1 + weight2*Im2); xlabel('Weighted Addition');

% Image Subtraction
background1 = imopen(Im1,strel('disk',15));
ImSub1 = imsubtract(Im1,background1);

background2 = imopen(Im1,strel('diamond',15));
ImSub2 = imsubtract(Im1,background2);

background3 = imopen(Im1,strel('octagon',15));
ImSub3 = imsubtract(Im1,background3);

background4 = imopen(Im1,strel('square',15));
ImSub4 = imsubtract(Im1,background4);

background5 = imopen(Im1,strel('cube',15));
ImSub5 = imsubtract(Im1,background5);

background6 = imopen(Im1,strel('sphere',15));
ImSub6 = imsubtract(Im1,background6);

figure(5)
subplot(2, 4, 1); imshow(Im1); xlabel('Original');
subplot(2, 4, 2); imshow(ImSub1); xlabel('Disk'); title('Image Subtraction');
subplot(2, 4, 3); imshow(ImSub2); xlabel('Diamond');
subplot(2, 4, 4); imshow(ImSub3); xlabel('Octogon');
subplot(2, 4, 5); imshow(ImSub4); xlabel('Square');
subplot(2, 4, 6); imshow(ImSub5); xlabel('Cube');
subplot(2, 4, 7); imshow(ImSub6); xlabel('Sphere');


% Add bias to image (Addition and Subtraction)
bright1 = 75; 
bright2 = 50; 

figure(6) 
I1_addBias = Im1 + bright1;
I1_subBias = Im1 - bright1;

I2_addBias = Im1 + bright2;
I2_subBias = Im1 - bright2;

subplot(2,3,2); imshow(Im1); xlabel('Original Greyscale Image, I');
title('Brightness Change');
subplot(2,3,1); imshow(I1_addBias); xlabel('I + bias(75)');
subplot(2,3,3); imshow(I1_subBias); xlabel('I - bias(75)');

subplot(2,3,5); imshow(Im1); xlabel('Original Greyscale Image, I');
subplot(2,3,4); imshow(I2_addBias); xlabel('I + bias(50)');
subplot(2,3,6); imshow(I2_subBias); xlabel('I - bias(50)');


% Image Multiplication & Division
alpha = 1.2; 
beta = 2;
I1_multiplyAlpha = alpha*Im1;
I1_multiplyBeta = Im1/beta;

figure(7)
subplot(2,3,2); imshow(Im1); xlabel('Original Greyscale Image, I');
title('Contrast Change');
subplot(2,3,1); imshow(I1_multiplyAlpha); xlabel('\alpha * I');
subplot(2,3,3); imshow(I1_multiplyBeta); xlabel('I / \beta');

% Image Transpose
figure(8)
subplot(1,3,1); imshow(Im1); xlabel('Original');
subplot(1,3,2); imshow(Im1(end:-1:1,:)); title('Image Transpose');
subplot(1,3,3); imshow(imrotate(flip(Im1,2),90));

% Binary
B1 = graythresh(Im1);
B2 = graythresh(Im2);
BW1 = imbinarize(Im1,B1);
BW2 = imbinarize(Im2,B2);

figure(9)
subplot(2, 2, 1); imshow(BW1); title('Binary Images');
subplot(2, 2, 2); imshow(BW2);

% Logical Operation
figure(10)
logAnd = bitand(Im1, Im2);
subplot(2, 3, 1); imshow(logAnd); xlabel('And')

logOr = bitor(Im1, Im2);
subplot(2, 3, 2); imshow(logOr); xlabel('Or'); title('Logical Operations');

logNot1 = bitcmp(Im1);
subplot(2, 3, 3); imshow(logNot1); xlabel('Not (Im1)');

logNot2 = bitcmp(Im2);
subplot(2, 3, 4); imshow(logNot2); xlabel('Not (Im2)');

logXor = bitxor(Im1, Im2);
subplot(2, 3, 5); imshow(logXor); xlabel('Xor');
