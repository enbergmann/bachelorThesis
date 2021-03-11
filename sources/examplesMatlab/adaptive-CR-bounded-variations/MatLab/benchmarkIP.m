function benchmarkIP(theta, minDof, alpha, gamma, epsilon, s)
    addpath(genpath(pwd));

    %% load image
    img = imread('cameraman.tif');
    imgSize = size(img);
%     % add 10 pixel frame of zeros
%     img = [zeros(imgSize(1),10), img, zeros(imgSize(1),10)];
%     img = [zeros(10,imgSize(2)+20); img; zeros(10,imgSize(2)+20)];
%     imgSize = imgSize + 20;

    % convert utf8 to double
    img = im2double(img);
    
    for j = 1:25
        img(j,:) = (j-1)*img(j,:)/25;
        img(imgSize(1)-(j-1),:) = (j-1)*img(imgSize(1)-(j-1),:)/25;
        img([j:imgSize(1)-(j-1)],j) = (j-1)*img([j:imgSize(1)-(j-1)],j)/25;
        img([j:imgSize(1)-(j-1)],imgSize(2)-(j-1)) ...
            = (j-1)*img([j:imgSize(1)-(j-1)],imgSize(2)-(j-1))/25;
    end
    
    if nargin < 1
        theta = 0.5;
    end
    if nargin < 2
        minDof = 1000;
    end
    if nargin < 3
        alpha = 1000;
    end
    if nargin < 4
        gamma = 1;
    end
    if nargin < 5
        epsilon = 1e-5;
    end
    if nargin < 6
        s = 0.5;
    end
    
    % rescale
    img = img*alpha;
    
    degreeF = 0;
    
    [c4n, n4e, n4sDb, ~] = loadGeometry('Square');
    
    folderName = ['results2/ImgProcessing/alpha', num2str(alpha), 'theta', num2str(ceil(theta*100)), ...
                  'gamma', num2str(ceil(gamma*100)), 'dof', num2str(minDof)];
    
    [meanUCR4e, c4n, n4e] = afemBV(@(x)f(x, img, imgSize), c4n, n4e, n4sDb, theta, minDof,...
                alpha, gamma, epsilon, s, degreeF, folderName);

    plotGreyScale(meanUCR4e, c4n, n4e, folderName);
end

function y = f(x, img, imgSize)
    nrPts = size(x,1);
    y = zeros(nrPts,1);
    
    for j = 1:nrPts
        curPt = x(j,:);
        pos1 = max(ceil(curPt(1)*imgSize(1)), 1);
        pos2 = max(ceil(curPt(2)*imgSize(2)), 1);
        y(j) = img(imgSize(2) - (pos2 - 1), pos1);
    end
end
