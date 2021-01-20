format long g;
format compact;
originalPic = imread('PandaOriginal(1).bmp'); %original picture
figure,imshow(originalPic);
title(sprintf('Image without noise (Original)'));
noisePic = imread('PandaNoise(1).bmp'); % picture with noise
[rows,columns, depth] = size(noisePic); %storing rows column and depth data from noise image
noiseError = immse(noisePic,originalPic);
noisy = imshow(noisePic);
title(sprintf('Noisy Image MSE: %0.5f', noiseError));
%apply fast fourier transform
fastfourier = fft2(double(noisePic));
frequencyImage = fftshift(fastfourier);
spectrum = log(1+abs(frequencyImage));
%plot spectrum to figure out threshold
plot(spectrum);
figure, imshow(spectrum,[]);
minVal = min(min(spectrum));
maxVal = max(max(spectrum));
%threshold set as 11 calculated from spectrum plot
threshold = 11;
%mask set to any value in spectrum which is above threshold
%notch filter applied in frequency domain
mask = spectrum > threshold;
%centre of spectrum excluded from being set to 0 to ensure main image isn'truined
mask(round((rows/2))-30:round((rows/2))+30, round((columns/2))-50:round((columns/2)+50)) = 0;
%bright spikes zeroed out
frequencyImage(mask) = 0;
figure,imshow(mask);
fftMasked = log(1+abs(frequencyImage));
minValue = min(min(fftMasked));
maxValue = max(max(fftMasked));
figure,imshow(fftMasked, [minValue maxValue]);
plot(log(abs(frequencyImage)));
%fft inversed back to spatial domain
inverseImage = ifft2(fftshift(frequencyImage));
filteredImage = abs(inverseImage);
minValue = min(min(filteredImage));
maxValue = max(max(filteredImage));
notchFilteredImage = imshow(filteredImage, [minValue maxValue]);
fftError = immse(uint8(getimage(notchFilteredImage)),originalPic);
title(sprintf('Notch Filtered Image MSE: %0.5f', fftError));
%low pass filter design
circle=55; % cut off frequency (circle) size
[m, n]=size(noisePic); %matrix storing the size of the noisey picture
spikes=zeros(m,n); %any spikes will be set to 0 to remove noise
for x=1:m
for y=1:n
if (x-m/2)^2 + (y-n/2)^2 < circle^2 % if number in matrix isnt in circle size
spikes(x,y) = 1; %all of the points with low frequency are set to 1
end
end
end
figure,imshow(spikes, []);

spectrum = frequencyImage; % spectrum(x,y) * spikes(x,y) to give newSpectrum(x,y)
newSpectrum = spectrum.*spikes; % convolution matrix
lowPassFilter = abs(ifft2(fftshift(newSpectrum))); % inverse back to spatial
figure,lowPassFilteredImage = imshow(lowPassFilter, []);
lowPassFilterError = immse(uint8(getimage(lowPassFilteredImage)), originalPic);
title(sprintf('low pass filter Image MSE: %0.5f', lowPassFilterError));
% mean filter (average) on cropped adaptive noise filtered image
cropOriginal = imcrop(uint8(originalPic),[140 40 340 316]);
figure,croppedOriginalImage = imshow(cropOriginal);
croppedOriginalHandle = image(cropOriginal);
title('BF-1 Cropped Original Image');
% cropping the Notch filtered image
cropFiltered = imcrop(uint8(getimage(notchFilteredImage)),[140 40 340 316]);
figure,cropFilteredImage = imshow(cropFiltered);
% mean squared error between the cropped filtered, and original images
test = uint8(getimage(cropFilteredImage));
test2 = croppedOriginalImage;
cropError = immse(uint8(getimage(cropFilteredImage)), cropOriginal);
title(sprintf('Cropped NF Image MSE: %0.3f', cropError));
% applying the mean filter
meanFilter = imboxfilt(uint8(getimage(cropFilteredImage)), 3); % 3x3 kernel
figure,bFilteredImage = imshow(meanFilter);
% mean squared error between the cropped mean filtered, and original images
meanFilterError = immse(uint8(getimage(bFilteredImage)), cropOriginal);
title(sprintf('Cropped M&NF MSE: %0.3f', meanFilterError));
%applying median filter in on notch filtered image spatial domain
median = medfilt2(uint8(getimage(notchFilteredImage)));
figure,medFilteredImage = imshow(median);
medianError = immse(uint8(getimage(medFilteredImage)), originalPic);
title(sprintf('Median and Notch Filtered Image MSE: %0.5f', medianError));
%applying gaussian noise removal on notch filtered image in spatial domain
gaussianFilter = imgaussfilt(uint8(getimage(medFilteredImage)),2);
figure,gaussianFilteredImage = imshow(gaussianFilter);
gaussianFilterError = immse(uint8(getimage(gaussianFilteredImage)), originalPic);
title(sprintf('G and MED Image MSE: %0.5f', gaussianFilterError));
%applying adaptive noise removal on notch filtered image in spatial domain
adaptiveFilter = wiener2(uint8(getimage(medFilteredImage)));
figure,adaptiveFilteredImage = imshow(adaptiveFilter);
adaptiveFilterError = immse(uint8(getimage(adaptiveFilteredImage)), originalPic);
title(sprintf('AD & MED Image MSE: %0.5f', adaptiveFilterError));