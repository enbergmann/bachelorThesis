function val = rhsImg %rhsImg(x, img)
  img = imread('cameraman.tif');
  img = [zeros(200), 0.33*ones(200); 0.66*ones(200), ones(200)];
  imgSize = size(img);
%   % add 10 pixel frame of zeros
%   img = [zeros(imgSize(1),10), img, zeros(imgSize(1),10)];
%   img = [zeros(10,imgSize(2)+20); img; zeros(10,imgSize(2)+20)];
%   imgSize = imgSize + 20;

  % convert utf8 to double
  img = im2double(img);
  
  % add fade to black on the edges
  for j = 1:25
      img(j, :) = (j-1)*img(j, :)/25;
        % first 25 rows
      img(imgSize(1)-(j-1), :) = (j-1)*img(imgSize(1)-(j-1), :)/25;
        % last 25 rows
      img([j:imgSize(1)-(j-1)], j) = (j-1)*img([j:imgSize(1)-(j-1)], j)/25;
        % first 25 columns except already done first and last j rows
      img([j:imgSize(1)-(j-1)], imgSize(2)-(j-1)) ...
          = (j-1)*img([j:imgSize(1)-(j-1)], imgSize(2)-(j-1))/25;
        % last 25 columns except already done first and lastj rows
      % always divide by 25 since j is 25 at most
  end
  
  alpha = 1000; % same as in the algorithm
  
  % rescale (since bartels formulates the energy slightly differently)
  %img = img*alpha;
  
  addpath(genpath('../'))
  [c4n, n4e, n4sDb, ~] = loadGeometry('Square',2);

  % [meanUCR4e, c4n, n4e] = afemBV(@(x)f(x, img, imgSize), c4n, n4e, n4sDb, theta, minDof,...
  %             alpha, gamma, epsilon, s, degreeF, folderName);
  % plotGreyScale(meanUCR4e, c4n, n4e, folderName);
  keyboard
end

function y = f(x, img, imgSize)
  % prob only works on Square domain right now
  nrPts = size(x, 1);
  y = zeros(nrPts, 1);
  
  for j = 1:nrPts
    curPt = x(j, :);
    pos1 = max(ceil(curPt(1)*imgSize(1)), 1);
    pos2 = max(ceil(curPt(2)*imgSize(2)), 1);
    y(j) = img(imgSize(2) - (pos2 - 1), pos1);
  end
end
