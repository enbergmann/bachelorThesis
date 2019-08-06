function benchmarkExact(theta, minDof, geometry, alpha, epsilon, tau)
    addpath(genpath(pwd));
    
    if nargin < 1
        theta = 0.5;
    end
    if nargin < 2
        minDof = 500;
    end
    if nargin < 3
        geometry = 'BigSquare';
    end
    if nargin < 4
        alpha = 1;
    end
    if nargin < 5
        epsilon = 1e-6;
    end
    if nargin < 6
        tau = 0.5;
    end
    degreeF = 5;
    gamma = 1;
    delta = 1;
    f = @(x)g(x, alpha, delta);
    
    [c4n, n4e, n4sDb, ~] = loadGeometry(geometry);
    
    afemBV(f, c4n, n4e, n4sDb, theta, minDof, alpha, gamma, epsilon, tau,...
            degreeF, geometry, @(x)uExact(x,delta), @(x)DuExact(x,delta));
end

function val = g(x,alpha,delta)
    r = sqrt(sum(x.^2,2));
    ind1 = (0 <= r & r <= 1/6);
    ind2 = (1/6 <= r & r <= 1/3);
    ind3 = (1/3 <= r & r <= 1/2);
    ind4 = (1/2 <= r & r <= 5/6);
    ind5 = (5/6 <= r & r <= 1);
    
    val = zeros(size(x,1),1);
    val(ind1) = alpha - 12*(2 - 9*r(ind1));
    val(ind2) = alpha*(1 + (6*r(ind2) - 1).^delta) - 1./r(ind2);
    val(ind3) = 2*alpha + 6*pi*sin(pi*(6*r(ind3) - 2)) - cos(pi*(6*r(ind3) - 2))./r(ind3);
    val(ind4) = 2*alpha*(5/2 - 3*r(ind4)).^delta + 1./r(ind4);
    val(ind5) = -3*pi*sin(pi*(6*r(ind5) - 5)) + (1 + cos(pi*(6*r(ind5) - 5)))./(2*r(ind5));
end

function val = uExact(x, delta)
    r = sqrt(sum(x.^2,2));
    ind1 = (r < 1/6);
    ind2 = (1/6 <= r & r < 1/3);
    ind3 = (1/3 <= r & r < 1/2);
    ind4 = (1/2 <= r & r < 5/6);
    ind5 = (5/6 <= r & r <= 1);
    
    val = zeros(size(x,1),1);
    val(ind1) = 1;
    val(ind2) = 1 + (6*r(ind2) - 1).^delta;
    val(ind3) = 2;
    val(ind4) = 2*(5/2 - 3*r(ind4)).^delta;
    val(ind5) = 0;
end

% compute exact gradient of u
function val = DuExact(x, delta)
    r = sqrt(sum(x.^2,2));
    ind2 = (1/6 <= r & r < 1/3);
    ind4 = (1/2 <= r & r < 5/6);
    
    val = zeros(size(x,1),2);
    tmp = 6*delta*(6*r(ind2) - 1).^(delta - 1)./r(ind2);
    val(ind2,:) = [tmp.*x(ind2,1), tmp.*x(ind2,2)];
    
    tmp = -6*delta*(5/2 - 3*r(ind4)).^(delta - 1)./r(ind4);
    val(ind4,:) = [tmp.*x(ind4,1), tmp.*x(ind4,2)];
end