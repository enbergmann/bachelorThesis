function addBoundary2image(imageName, boundaryType)
  if nargin<2, boundaryType = 'f2b'; end;

  % read image and convert utf8 to double
  img = im2double(imread(imageName));
  imgSize = size(img);

  if strcmp(boundaryType, 'f2b') 
    % add fade to black on the edges for 25 pixels (0 boundary conditions)
    %
    % img(1, 1),     [0, 1] top left pixel
    % img(1, end),   [1, 1] top right pixel
    % img(end, 1),   [0, 0] bottom left pixel
    % img(end, end), [1, 0] bottom right pixel
    for j = 1:25
      img(j, j:imgSize(2)-(j-1)) = (j-1)*img(j, j:imgSize(2)-(j-1))/25;
        % first 25 rows above (including) main diagonals
      img(imgSize(1)-(j-1), j:imgSize(2)-(j-1)) = ...
        (j-1)*img(imgSize(1)-(j-1), j:imgSize(2)-(j-1))/25;
        % last 25 rows below (including) main diagonals
      img(j+1:imgSize(1)-j, j) = (j-1)*img(j+1:imgSize(1)-j, j)/25;
        % first 25 columns left of (excluding) main diagonals
      img(j+1:imgSize(1)-j, imgSize(2)-(j-1)) = ...
        (j-1)*img(j+1:imgSize(1)-j, imgSize(2)-(j-1))/25;
        % last 25 columns right of (excluding) main diagonals
    end
  elseif strcmp(boundaryType, 'fullb') 
    % add 10 pixel frame of zeros
    img = [zeros(imgSize(1),10), img, zeros(imgSize(1),10)];
    img = [zeros(10,imgSize(2)+20); img; zeros(10,imgSize(2)+20)];
  else
    error('boundaryType unknown');
  end

  imwrite(img, sprintf('%s%s', boundaryType, imageName));
end
