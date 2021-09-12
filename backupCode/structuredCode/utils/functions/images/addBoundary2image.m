function addBoundary2image(imageName, boundaryType)
%% DOC
% Adds a boundary to a given image guaranteeing that the value of the outermost
% pixels is 0 (black for a grayscale picture).
%
% addBoundary2image.m
% input: imageName    - 'char array with exactly one row' containing the path
%                       (including filename and file extension) of the given
%                       image
%        boundaryType - 'char array with exactly one row' containing the type
%                       of boundary that must be added to the picture
%                       (optional, default value is 'f2b')
%          Options:
%              'f2b': fade to black towards the edge during the outermost 25
%                     pixels of the image
%            'fullb': add 10 new outermost pixels to the image that are
%                     completely black

%% INIT
  if nargin<2, boundaryType = 'f2b'; end

  % read image and convert utf8 to double
  img = im2double(imread(imageName));
  [imgSizeX, imgSizeY] = size(img);

%% MAIN
  if strcmp(boundaryType, 'f2b') 
    % add fade to black on the edges for 25 pixels (0 boundary conditions)
    %
    % img(1, 1),     [0, 1] top left pixel
    % img(1, end),   [1, 1] top right pixel
    % img(end, 1),   [0, 0] bottom left pixel
    % img(end, end), [1, 0] bottom right pixel
    for j = 1:25
      img(j, j:imgSizeY - (j - 1)) = ...
        (j - 1)*img(j, j:imgSizeY - (j - 1))/25;
        % first 25 rows above (including) main diagonals
      img(imgSizeX - (j - 1), j:imgSizeY - (j - 1)) = ...
        (j - 1)*img(imgSizeX - (j - 1), j:imgSizeY - (j - 1))/25;
        % last 25 rows below (including) main diagonals
      img(j + 1:imgSizeX - j, j) = (j - 1)*img(j + 1:imgSizeX - j, j)/25;
        % first 25 columns left of (excluding) main diagonals
      img(j + 1:imgSizeX - j, imgSizeY - (j - 1)) = ...
        (j - 1)*img(j + 1:imgSizeX - j, imgSizeY - (j - 1))/25;
        % last 25 columns right of (excluding) main diagonals
    end
  elseif strcmp(boundaryType, 'fullb') 
    % add 10 pixel frame of zeros
    img = [zeros(imgSizeX, 10), img, zeros(imgSizeX, 10)];
    img = [zeros(10, imgSizeY + 20); img; zeros(10, imgSizeY + 20)];
  else
    error('boundaryType unknown');
  end

  % save resulting image
  imwrite(img, sprintf('%s%s', boundaryType, imageName));
end
