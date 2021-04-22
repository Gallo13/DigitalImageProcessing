% Created by Jessica Gallo
% Date CreatedL 3/22/2021
% Last Modified: 3/22/2021
% CSC741 Digital Image Processing
% Prof. Agaian
% Lecture 7 HW

% To do:
% - get total PSRN for each bitplane
% - fix noise estimation from last time

close all; clear; clc;

% Load Image, Greyscale, Resize
%Im1 = imread('thermalImage3.jpg');
Im1 = imread('thermalImage3.jpg');
Im1G = rgb2gray(Im1);
%Im1G = imresize(120, 120);

% Bit Plane Slicing Pogram
% Apply bit-plane slicing in Matlab to read cameraman image, then extract 
% the image of bit 6
% ========================================================================
Im1G_double = double(Im1G);

b1 = mod(Im1G_double, 2);
b2 = mod(floor(Im1G_double/2), 2);
b3 = mod(floor(Im1G_double/4), 2);
b4 = mod(floor(Im1G_double/8), 2);
b5 = mod(floor(Im1G_double/16), 2);
b6 = mod(floor(Im1G_double/32), 2);
b7 = mod(floor(Im1G_double/64), 2);
b8 = mod(floor(Im1G_double/128), 2);

Im1G_combined = (2 * (2 * (2 * (2 * (2 * (2 * (2 * b8 + b7) + b6) + b5) + b4) + b3) + b2) + b1);
% Im1G_combined = (2 * (2 * (2 *  b3) + b2) + b1);

figure(1);

subplot(2, 5, 1); imshow(Im1G); xlabel('Original');
subplot(2, 5, 2); imshow(b1); xlabel('Bit Plane 1');
subplot(2, 5, 3); imshow(b2); xlabel('Bit Plane 2');
subplot(2, 5, 4); imshow(b3); xlabel('Bit Plane 3');
subplot(2, 5, 5); imshow(b4); xlabel('Bit Plane 4');
subplot(2, 5, 6); imshow(b5); xlabel('Bit Plane 5');
subplot(2, 5, 7); imshow(b6); xlabel('Bit Plane 6');
subplot(2, 5, 8); imshow(b7); xlabel('Bit Plane 7');
subplot(2, 5, 9); imshow(b8); xlabel('Bit Plane 8');
subplot(2, 5, 10); imshow(uint8(Im1G_combined)); xlabel('Recombined Image');

% subplot(2, 1, 1); imshow(uint8(Im1G_combined)); xlabel('Recombined Image 3-8');

% PSNR (Peak Signal to Noise Ratio)
Im1G_combined = uint8(Im1G_combined);
Im1G = uint8(Im1G);
[peaksnr0, snr0] = psnr(Im1G_combined, Im1G);
fprintf('\nThe PSNR between the original image and the combined bit plane image is %0.4f', peaksnr0);
fprintf('\nThe SNR between the original image and the combined bit plane image is %0.4f', snr0);

% Develop a tool to decompose an image into bit plane. Use the AND 
% operator between slices and display your results.
% ========================================================================
Im1G = uint8(Im1G);
[row, col] = size(Im1G);
for i=1:1:row
    for j=1:1:col
        MSB(i,j)=bitand(Im1G(i,j), bin2dec('10000000'));
        LSB(i,j)=bitand(Im1G(i,j), bin2dec('00000001'));
        b22(i,j)=bitand(Im1G(i,j), bin2dec('01000000'));
        b33(i,j)=bitand(Im1G(i,j), bin2dec('00100000'));
        b44(i,j)=bitand(Im1G(i,j), bin2dec('00010000'));
        b55(i,j)=bitand(Im1G(i,j), bin2dec('00001000'));
        b66(i,j)=bitand(Im1G(i,j), bin2dec('00000100'));
        b77(i,j)=bitand(Im1G(i,j), bin2dec('00000010'));
    end
end

Im1G_combined2 = (2 * (2 * (2 * (2 * (2 * (2 * (2 * MSB + b77) + b66) + b55) + b44) + b33) + b22) + LSB);
%Im1G_combined = (2 * (2 * (2 * MSB) + b3) + b2);

figure(2);
subplot(2, 4, 1); imshow(MSB); xlabel('MSB');
subplot(2, 4, 2); imshow(LSB); xlabel('LSB');
subplot(2, 4, 3); imshow(b22); xlabel('Bit plane 2');
subplot(2, 4, 4); imshow(b33); xlabel('Bit plane 3');
subplot(2, 4, 5); imshow(b44); xlabel('Bit plane 4');
subplot(2, 4, 6); imshow(b55); xlabel('Bit plane 5');
subplot(2, 4, 7); imshow(b66); xlabel('Bit plane 6');
subplot(2, 4, 8); imshow(b77); xlabel('Bit plane 7');
%figure(3);
%subplot(2, 1, 1); imshow(Im1G_combined2);



% Develop a tool to filter noisy images using bit plane decomposition
% ========================================================================

Im1G_noise1 = imnoise(Im1G, 'salt & pepper', 0.02); % gaussian, 0 mean, 0.02 variance

Im1G_noise1_double = double(Im1G_noise1);

b1n = mod(Im1G_noise1_double, 2);
b1n = medfilt2(b1n);
b2n = mod(floor(Im1G_noise1_double/2), 2);
b2n = medfilt2(b2n);
b3n = mod(floor(Im1G_noise1_double/4), 2);
b3n = medfilt2(b3n);
b4n = mod(floor(Im1G_noise1_double/8), 2);
b4n = medfilt2(b4n);
b5n = mod(floor(Im1G_noise1_double/16), 2);
b5n = medfilt2(b5n);
b6n = mod(floor(Im1G_noise1_double/32), 2);
b6n = medfilt2(b6n);
b7n = mod(floor(Im1G_noise1_double/64), 2);
b7n = medfilt2(b7n);
b8n = mod(floor(Im1G_noise1_double/128), 2);
b8n = medfilt2(b8n);

Im1G_noise1_combined = (2 * (2 * (2 * (2 * (2 * (2 * (2 * b8n + b7n) + b6n) + b5n) + b4n) + b3n) + b2n) + b1n);
% Im1G_combined = (2 * (2 * (2 *  b3) + b2) + b1);

figure(1);

subplot(2, 5, 1); imshow(Im1G_noise1); xlabel('Original Noise');
subplot(2, 5, 2); imshow(b1n); xlabel('Bit Plane 1');
subplot(2, 5, 3); imshow(b2n); xlabel('Bit Plane 2');
subplot(2, 5, 4); imshow(b3n); xlabel('Bit Plane 3');
subplot(2, 5, 5); imshow(b4n); xlabel('Bit Plane 4');
subplot(2, 5, 6); imshow(b5n); xlabel('Bit Plane 5');
subplot(2, 5, 7); imshow(b6n); xlabel('Bit Plane 6');
subplot(2, 5, 8); imshow(b7n); xlabel('Bit Plane 7');
subplot(2, 5, 9); imshow(b8n); xlabel('Bit Plane 8');
subplot(2, 5, 10); imshow(uint8(Im1G_noise1_combined)); xlabel('Recombined Image');

% subplot(2, 1, 1); imshow(uint8(Im1G_combined)); xlabel('Recombined Image 3-8');

% PSNR (Peak Signal to Noise Ratio)
Im1G_noise1_combined = uint8(Im1G_noise1_combined);
Im1G = uint8(Im1G);
[peaksnr0, snr0] = psnr(Im1G_noise1_combined, Im1G);
fprintf('\nThe PSNR between the original image and the noisy combined bit plane image is %0.4f', peaksnr0);
fprintf('\nThe SNR between the original image and the noisy combined bit plane image is %0.4f', snr0);

Im1G_noise1_combined = uint8(Im1G_noise1_combined);
Im1G = uint8(Im1G);
[peaksnr0, snr0] = psnr(Im1G_noise1_combined, Im1G_noise1);
fprintf('\nThe PSNR between the original noisy image and the noisy combined bit plane image is %0.4f', peaksnr0);
fprintf('\nThe SNR between the original noisy image and the noisy combined bit plane image is %0.4f', snr0);

% Develop a tool to encrypt images using bit plane decomposition.
% ========================================================================
% ref: https://www.geeksforgeeks.org/lsb-based-image-steganography-using-matlab/
% Encryption
% -------------------
% Message to be embedded
message='Hello World!'; 
% Length of the message where each character is 8 bits
len = length(message) * 8;
% Get all the ASCII values of the characters of the message
ascii_value = uint8(message); 
% Convert the decimal values to binary
bin_message = transpose(dec2bin(ascii_value, 8));
% Get all the binary digits in separate row
bin_message = bin_message(:);
% Length of the binary message
N = length(bin_message);
% Converting the char array to numeric array
bin_num_message=str2num(bin_message);
% Initialize output as input
output = Im1G;
  
% Get height and width for traversing through the image
height = size(Im1G, 1);
width = size(Im1G, 2);
  
% Counter for number of embedded bits
embed_counter = 1;
  
% Traverse through the image
for i = 1 : height
    for j = 1 : width
        % If more bits are remaining to embed
        if(embed_counter <= len)
            % Finding the Least Significant Bit of the current pixel
            LSB = mod(double(Im1G(i, j)), 2);
            % Find whether the bit is same or needs to change
            temp = double(xor(LSB, bin_num_message(embed_counter)));
            % Updating the output to input + temp
            output(i, j) = Im1G(i, j)+temp;
            % Increment the embed counter
            embed_counter = embed_counter+1;
        end
    end
end
  
% Write both the input and output images to local storage
% Mention the path to a folder here.
figure(4);
subplot(2, 1, 1); imshow(Im1G); xlabel('Original');
subplot(2, 1, 2); imshow(output); xlabel('Encrypted Image');

% Decryption
% -------------------
% https://github.com/shamilee05/Steganography-LSB/blob/master/Steganography.m
count=1;
message_in_bits='';
for i=1:height
    for j=1:width
        %For all the characters in the message
        if count<=len
            
            %Retrieve the LSB of the intensity level of the pixel
            LSB=mod(output(i,j),2);
            
            %Append into message_in_bits to get bit sequence of message
            message_in_bits=append(message_in_bits,num2str(LSB));
            
            count=count+1;
        end
    end
end

%Converting the bit sequence into the original message
i=1;
original_message='';
while i<=len
    %Take a set of 8 bits at a time
    %Convert the set of bits to a decimal number
    %Convert the decimal number which is the ascii value to its corresponding character
    %Append the obtained character into the resultant string
    original_message=append(original_message,char(bin2dec(message_in_bits(1,i:i+7))));
    i=i+8;
end

disp(['The original message is: ',original_message]);

% In this task you need to load an image file into MATLAB and slice it 
% bitwise into 8 planes, then try to compress it by eliminating the 
% unnecessary bits to reduce the whole image size. Practically you need to 
% implement a code that compresses an image using bit plane decomposition.
% https://shubbakom.files.wordpress.com/2014/03/imageprocessing-project3.pdf
% ========================================================================

% size
Im1G = uint8(Im1G);
[row, col] = size(Im1G);
% bitwise seperation
Inew = zeros(row, col);
mask = 1;
for a = 1:8
	for i = 1:1:row
		for j = 1:1:col
			Inew(i, j) = bitand(Im1G(i,j),mask);
		end
	end
	mask = mask*2;
	figure, imshow(Inew);
end

%Compressing image
% by eliminating the least significant four bits then the picture will 
% need only 4 bits (16 levels) instead of 8 bits (256 levels) ie half 
% storage size)
mask = 240;
for i = 1:1:row
	for j = 1:1:col
		Inew(i, j) = bitand(Im1G(i,j),mask);
	end
end
figure; imshow(uint8(Inew))

% the difference
Idif = Im1G - uint8(Inew);
figure; imagesc(Idif), colormap grey

% Estimate noise parameters from a single image
% ================================================================================================================
% Segment a small patch of image with constant grey level
Im1G_noise1 = imnoise(Im1G, 'gaussian', 0.02); % gaussian, 0 mean, 0.02 variance
% light in face
Im1G_noise1_GC = grayconnected(Im1G_noise1,50,50); % segments the image at coordinates 50,50
% dark around face
% Im1G_noise1_GC = grayconnected(Im1G_noise1,10,110); % segments the image at coordinates 50,50

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
% Compute the mean and variance of small patch of noisy image
x = imhist(Im1G_noise1_GC);
Im1G_noise1_GC_mean = mean2(x);
fprintf('The mean of the small patch of the noise image is %0.4d', Im1G_noise1_GC_mean);
Im1G_noise1_GC_var = var(x);
fprintf('\nThe variance of the small patch of the noise image is %0.4d', Im1G_noise1_GC_var);
Im1G_noise1_GC_std = std(x);
fprintf('\nThe standard deviation of the small patch of the noise image is %0.4d', Im1G_noise1_GC_std);
% Computer mean and variance of original image
x = imhist(Im1G);
Im1G_mean = mean2(x);
fprintf('\nThe mean of the original image is %0.4d', Im1G_mean);
Im1G_var = var(x);
fprintf('\nThe variance of the original image is %0.4d', Im1G_var);
Im1G_std = std(x);
fprintf('\nThe standard deviation of the original image is %0.4d', Im1G_std);
% Computer mean and variance of original noisy image
x = imhist(Im1G_noise1);
Im1G_noise1_mean = mean2(x);
fprintf('\nThe mean of the original noise image is %0.4d', Im1G_noise1_mean);
Im1G_noise1_var = var(x);
fprintf('\nThe variance of the original noise image is %0.4d', Im1G_noise1_var);
Im1G_noise1_std = std(x);
fprintf('\nThe standard deviation of the original noise image is %0.4d', Im1G_noise1_std);

