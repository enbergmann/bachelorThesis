function y = rhsImg(x, img, imgSize)
% TODO after this works create denoise exampe, see if the algorithm denoises
% TODO benchmark needs modes like (function, image, denoise)
% TODO comment this when it's beautiful
% TODO name is prob more sth like 'projectImageToSquare'
  % prob only works on Square domain right now
  
  % img(1, 1),     [0, 1] top left pixel
  % img(1, end),   [1, 1] top right pixel
  % img(end, 1),   [0, 0] bottom left pixel
  % img(end, end), [1, 0] bottom right pixel
  
  y = img(sub2ind(imgSize, ...
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
