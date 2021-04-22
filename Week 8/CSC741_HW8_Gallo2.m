% Created by: Jessica Gallo
% Date Created: 4/5/2021
% Last Modified: 4/8/2021
% HW 8
% Fourier Transforms

% Encrypt
% =========================================================================
% ref: https://www.et.byu.edu/~vps/ME505/AAEM/V5-01.pdf

clc; % Clear the command window.
close all; % Close all figures (except those of imtool.)
clear; % Erase all existing variables. Or clearvars if you want.

%2D FFT Demo
%Import image
image = im2double(imread('00483_s_20aqapbvgk0483.png','png'));
image = image(:,:,1);
image = rescale(image);

%Display image
figure(1)
imshow(image)

%Perform 2D FFTs
imageFFT = fft2(image);
imageFFTMag = abs(imageFFT);
imgaeFFTAng = angle(imageFFT);

%Display magnitude and phase of 2D FFTs
figure(2)
imshow(log(1 + abs(fftshift(imageFFT))), []);
colormap gray

%Save the imageFFTMag to verify results at end
diff1 = imageFFTMag;

%Encrypt by randomization
rows = zeros(length(imageFFTMag(1,:)),1);
for index = 1:length(imageFFTMag(1,:))
 rows(index) = randi([0 2*1000])-1000;
end
cols = zeros(length(imageFFTMag(:,1)),1);
for index = 1:length(imageFFTMag(:,1))
 cols(index) = randi([0 2*1000])-1000;
end
for index = 1:length(rows)
 imageFFTMag(:,index) = circshift(imageFFTMag(:,index),[rows(index) 0]);
end
for index = 1:length(cols)
 imageFFTMag(index,:) = circshift(imageFFTMag(index,:),[0 cols(index)]);
end

%Create new image
imageCryptFFT = imageFFTMag.*exp(1i*imgaeFFTAng);

%Display cyrpt magnitude and phase of 2D FFTs
figure(3)
imshow(log(1 + abs(fftshift(imageCryptFFT))), []);

%Perform inverse 2D FFTs on image
imageCrypt = ifft2(imageCryptFFT);
imageCrypt = [real(imageCrypt), imag(imageCrypt)];

%Save the key. This will allow the image to be decrypted
key = [rows; cols];
fileID = fopen('key2.txt','w');
fprintf(fileID,'%d\n',key);
fclose(fileID);

%We have to separate real and imaginary, because imwrite will only write real values
imageCryptWrite = rescale(im2double(imageCrypt));
imageCryptWriteU8 = im2uint8(imageCryptWrite);

%Save Crypt image
imwrite(imageCryptWriteU8,'imagecrypt2.png','png');

% Dencrypt
% =========================================================================
% ref: https://www.et.byu.edu/~vps/ME505/AAEM/V5-01.pdf

%Import image
imageCrypt = im2double(imread('imagecrypt2.png','png'));

%read in key
fileID = fopen('key2.txt','r');
key = fscanf(fileID,'%f');

%Display image
figure(4)
imshow(imageCrypt)

%Join the complex and real parts
imageCryptreal = imageCrypt(:,1:length(imageCrypt(1,:))/2);
imageCryptimag = imageCrypt(:,length(imageCrypt(1,:))/2+1:end);
imageCrypt = complex(imageCryptreal,imageCryptimag);

%Display image
figure(5)
imshow(rescale(im2double(abs(imageCrypt))))

%Perform 2D FFTs
imageCryptFFT = fft2(imageCrypt);
imageCryptFFTMag = abs(imageCryptFFT);
imageCryptFFTAng = angle(imageCryptFFT);

%Display magnitude and phase of 2D FFTs
figure(6)
imshow(log(1 + abs(fftshift(imageCryptFFT))), []);

%Decrypt
keyindex = 1;
key = flipud(key);
for index = 1:length(imageCryptFFTMag(:,1))
 imageCryptFFTMag(length(imageCryptFFTMag(:,1))-index+1,:) = circshift(imageCryptFFTMag(length(imageCryptFFTMag(:,1))-index+1,:),[0 -key(keyindex)]);
 keyindex = keyindex + 1;
end
for index = 1:length(imageCryptFFTMag(1,:))
 imageCryptFFTMag(:,length(imageCryptFFTMag(1,:))-index+1) = circshift(imageCryptFFTMag(:,length(imageCryptFFTMag(1,:))-index+1),[-key(keyindex) 0]);
 keyindex = keyindex + 1;
end

%Save the imageFFTMag to verify results at end
diff2 = imageCryptFFTMag;

%Create new image
imageFFT = imageCryptFFTMag.*exp(1i*imageCryptFFTAng);

%Display cyrpt magnitude and phase of 2D FFTs
figure(7)
imshow(log(1 + abs(fftshift(imageFFT))), []);

%Perform inverse 2D FFTs on image
image = ifft2(imageFFT);

%Calculate limits for plotting
imageMin = min(min(min(abs(image))));
imageMax = max(max(max(abs(image))));

%Display cyrpt images
figure(8)
subplot(2, 1, 1); imshow(abs(image), [imageMin imageMax]);
subplot(2, 1, 2); imshow(image)
colormap gray

%Verify results
diff = abs(diff1) - abs(diff2);
figure(9)
plot(diff)

%diff3 = abs(image) - abs(image);

% add mse and psnr scores

% MSE
% original thermal image & decrypted thermal image
err1 = immse(image, imageCrypt);
fprintf('\n The MSE of the thermal image is %0.4f\n', err1);


err0 = immse(diff1, diff2);
fprintf('\n The MSE of the thermal image is %0.4f\n', err0);


% original car image & decrypted car image
err2 = immse(image, imageCrypt);
fprintf('\n The MSE of the car image is %0.4f\n', err2);

% PSNR (Peak Signal to Noise Ratio)
[peaksnr0, snr0] = psnr(image, imageCrypt);
fprintf('\nThe PSNR between the original image and the decrypted image is %0.4f', peaksnr0);
fprintf('\nThe SNR between the original image and the decrypted image is %0.4f', snr0);

[peaksnr01, snr01] = psnr(diff1, diff2);
fprintf('\nThe PSNR between the original image and the decrypted image is %0.4f', peaksnr01);
fprintf('\nThe SNR between the original image and the decrypted image is %0.4f', snr01);

% Filtering
% ========================================================================

filteredImage = medfilt2(image);
filteredImage2 = medfilt2(imageCrypt);

figure(10)
subplot(2, 1, 1); imshow(image);
subplot(2, 1, 2); imshow(imageCrypt);

%{
clc; 
clear;
close all;
% Filtering
% =========================================================================
I = imread('thermalImage3.jpg');
I=I-mean(I(:));
f = fftshift(fft2(I));
fabs=abs(f);
figure(1)
subplot(1,1,1)
imshow(fabs,[])
thresh=0.9*max(fabs(:));
% Find pixels that are brighter than the threshold.
mask = fabs > thresh; 
% Erase those from the image
fabs(mask) = 0;
% Shift back and inverse fft
filteredImage = ifft2(fftshift(fabs)) + mean2(I);
figure(2)

imshow(filteredImage, []);
%}