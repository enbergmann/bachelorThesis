function f = image2function(imageName, parAlpha)
%% DOC
% Produces an input signal for the nonconforming problem on the unit square
% from a given image.
%
% image2function.m
% input: imageName - 'char array with exactly one row' containing the path
%                    (including filename and file extension) of the given image
%        parAlpha  - 'double' containing the parmeter \alpha from the problem
%
% output: f - 'function_handle' of the input signal

%% INIT
  % read image and convert utf8 to double
  img = im2double(imread(imageName));
  
  % rescale (since f = \alpha g for input signal f of the nonconforming problem
  % and input signal g of the ROF model problem)
  img = img*parAlpha;

  % convert image to function
  % img(1, 1),     [0, 1] top left pixel
  % img(1, end),   [1, 1] top right pixel
  % img(end, 1),   [0, 0] bottom left pixel
  % img(end, end), [1, 0] bottom right pixel
  imgSize = size(img);

%% MAIN
  f = @(x) img(sub2ind(imgSize, ...
    max(round(imgSize(1) - x(:, 2)*imgSize(1)), 1), ...
    max(round(x(:, 1)*imgSize(2)), 1)));
end
