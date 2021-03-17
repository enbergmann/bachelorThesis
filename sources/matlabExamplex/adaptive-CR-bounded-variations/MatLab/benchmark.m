function benchmark(theta, minDof, geometry, alpha, epsilon, tau)
    addpath(genpath(pwd));
    
    if nargin < 1
        theta = 0.5;
    end
    if nargin < 2
        minDof = 1000;
    end
    if nargin < 3
        geometry = 'BigSquare';
    end
    if nargin < 4
        alpha = 1;
    end
    if nargin < 5
        epsilon = 1e-4;
    end
    if nargin < 6
        tau = 0.5;
    end
    
    degreeF = 0;
    gamma = 1;
    
    [c4n, n4e, n4sDb, ~] = loadGeometry(geometry);
    
    afemBV(@f, c4n, n4e, n4sDb, theta, minDof,...
                alpha, gamma, epsilon, tau, degreeF, geometry);
end

% function y = f(x)
%     y = ones(size(x,1), 1);
% end

function y = f(x)
    y = zeros(size(x,1), 1);
    ind = (-1/2 < x(:,1) & x(:,1) < 1/2 & -1/2 < x(:,2) & x(:,2) < 1/2);
    y(ind) = 100;
end