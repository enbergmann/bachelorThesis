function addNoise2image(imageName, snr)

  if nargin<2, snr = 20; end;

  % read image and convert utf8 to double
  img = im2double(imread(imageName));

  img = awgn(img, snr, 'measured');

  imwrite(img, sprintf('awgnSnr%d%s', snr, imageName));
end
