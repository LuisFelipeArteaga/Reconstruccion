clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 18;

button = menu('Use which demo image?', 'Coins', 'Moon', 'Circles', 'Spine');
if button == 1
	baseFileName = 'coins.png';
	thresholdValue = 80;
elseif button == 2
	baseFileName = 'moon.tif';
	thresholdValue = 80;
elseif button == 3
	baseFileName = 'circles.png';
	thresholdValue = 0.5;
elseif button == 4
	baseFileName = 'spine.tif';
	thresholdValue = 20;
end

% Read in a standard MATLAB gray scale demo image.
folder = fileparts(which('cameraman.tif')); % Determine where demo folder is (works with all versions).
% baseFileName = 'eight.tif';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
	% File doesn't exist -- didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
grayImage = imread(fullFileName);
% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
if numberOfColorBands > 1
	% It's not really gray scale like we expected - it's color.
	% Convert it to gray scale by taking only the green channel.
	grayImage = grayImage(:, :, 2); % Take green channel.
end
% Display the original gray scale image.
subplot(2, 3, 1);
imshow(grayImage, []);
axis on;
title('Original Grayscale Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
drawnow;

% Let's compute and display the histogram.
[pixelCount, grayLevels] = imhist(grayImage);
subplot(2, 3, 2); 
bar(grayLevels, pixelCount);
grid on;
title('Histogram of original image', 'FontSize', fontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.

% Construct an initial mask that is the convex hull of everything.
mask = grayImage > thresholdValue; % Use for spine.tif
% Display the image.
subplot(2, 3, 3);
imshow(mask, []);
axis on;
caption = sprintf('Thresholded at %.1f to get binary image', thresholdValue);
title(caption, 'FontSize', fontSize);

% Get an initial mask that is approximate and smoothed - a fair "first guess."
% Get the convex hull to get a smoother starting mask.
mask = bwconvhull(mask, 'Union');
% Alternate way to get a smoother mask using morphology.
% mask = imclose(mask, true(45)); % Use imclose, or imdilate but not both.
% mask = imdilate(mask, true(35));
% mask = imfill(mask, 'holes');
% Display the smoothed initial mask image.
subplot(2, 3, 4);
imshow(mask);
axis on;
title('Initial mask using convex hull', 'FontSize', fontSize);
drawnow;

% Now find the improved outer boundary using active contours.
numberOfIterations = 400;
bw = activecontour(grayImage, mask, numberOfIterations, 'edge');
subplot(2, 3, 5);
imshow(bw);
axis on;
caption = sprintf('Final outer boundary mask using %d iterations', numberOfIterations);
title(caption, 'FontSize', fontSize);
drawnow;

% Display the original gray scale image in the lower left
% So we can display boundaries over it.
subplot(2, 3, 6);
imshow(grayImage, []);
axis on;
hold on;

% Display the initial contour on the original image in blue.
% Display the final contour on the original image in red.
contour(mask,1,'b', 'LineWidth', 4); 
contour(bw, 1,'r', 'LineWidth', 4); 
title('Image with initial and final contours overlaid', 'FontSize', fontSize);

uiwait(msgbox('Done with activecontour demo'));
% close all;  % Close all figures (except those of imtool.)
