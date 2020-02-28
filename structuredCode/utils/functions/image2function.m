function f = image2function(imageName, parAlpha, addNoise)
  %TODO make this prettier at some point, but it works just fine for now
  % even though only on square
% TODO comment this when it's beautiful
% TODO name is prob more sth like 'projectImageToSquare'

  % prob only works on Square domain right now

  % TODO call function image2function or sth

  % read image and convert utf8 to double
  img = im2double(imread(imageName));
  if addNoise, img = imnoise(img, 'gaussian'); end
  imgSize = size(img);

  %   % add 10 pixel frame of zeros
  %   img = [zeros(imgSize(1),10), img, zeros(imgSize(1),10)];
  %   img = [zeros(10,imgSize(2)+20); img; zeros(10,imgSize(2)+20)];
  %   imgSize = imgSize + 20;
  
  % add fade to black on the edges for 25 pixels (0 boundary conditions)
  for j = 1:25
    img(j, :) = (j-1)*img(j, :)/25;
      % first 25 rows
    img(imgSize(1)-(j-1), :) = (j-1)*img(imgSize(1)-(j-1), :)/25;
      % last 25 rows
    img(j:imgSize(1)-(j-1), j) = (j-1)*img(j:imgSize(1)-(j-1), j)/25;
      % first 25 columns except already done first and last j rows
    img(j:imgSize(1)-(j-1), imgSize(2)-(j-1)) ...
        = (j-1)*img(j:imgSize(1)-(j-1), imgSize(2)-(j-1))/25;
      % last 25 columns except already done first and lastj rows
    % always divide by 25 since j is 25 at most
  end
    
  % rescale (since bartels formulates the energy slightly differently)
  img = img*parAlpha;

  % img(1, 1),     [0, 1] top left pixel
  % img(1, end),   [1, 1] top right pixel
  % img(end, 1),   [0, 0] bottom left pixel
  % img(end, end), [1, 0] bottom right pixel
  
  f = @(x) img(sub2ind(imgSize, ...
    max(round(imgSize(1) - x(:, 2)*imgSize(1)), 1), ...
    max(round(x(:, 1)*imgSize(2)), 1)));
  
  % nrPts = size(x, 1);
  % y = zeros(nrPts, 1);
  % 
  % for j = 1:nrPts
  %   curPt = x(j, :);
  %   y(j) = img(max(round(imgSize(1) - curPt(2)*imgSize(1)), 1), ...
  %     max(round(curPt(1)*imgSize(2)), 1));
  %   % y(j) = 1/4 * (...
  %   %   img(max(ceil(imgSize(1) - curPt(2)*imgSize(1)), 1), ...
  %   %   max(floor(curPt(1)*imgSize(2)), 1)) + ... %lower left pixel
  %   %   img(max(ceil(imgSize(1) - curPt(2)*imgSize(1)), 1), ...
  %   %   max(ceil(curPt(1)*imgSize(2)), 1)) + ... %lower right pixel
  %   %   img(max(floor(imgSize(1) - curPt(2)*imgSize(1)), 1), ...
  %   %   max(floor(curPt(1)*imgSize(2)), 1)) + ... %upper left pixel
  %   %   img(max(floor(imgSize(1) - curPt(2)*imgSize(1)), 1), ...
  %   %   max(ceil(curPt(1)*imgSize(2)), 1)));     %upper right pixel
  %   % NOTE this was the mean of all nearest pixel, one could furthermore
  %   % take into account how close the pixel are (postponed for now)
  %  end
end
