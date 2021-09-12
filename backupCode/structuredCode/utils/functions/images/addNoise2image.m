function addNoise2image(imageName, snr)
%% DOC
% Adds white Gaussian noise (AWGN) to an image.
%
% addNoise2image.m
% input: imageName - 'char array with exactly one row' containing the path
%                    (including filename and file extension) of the given image
%        snr       - 'double' containing signal-to-noise (SNR) ratio for
%                    MATLABs awgn function

%% INIT
  if nargin<2, snr = 20; end

  % read image and convert utf8 to double
  img = im2double(imread(imageName));

%% MAIN
  % add white Gaussian noise to image
  img = awgn(img, snr, 'measured');

  % save image
  imwrite(img, sprintf('awgnSnr%d%s', snr, imageName));
end
