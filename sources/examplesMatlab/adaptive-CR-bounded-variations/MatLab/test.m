function test
    addpath(genpath(pwd));
    
    [c4n, n4e, n4sDb, ~] = loadGeometry('Lshape', 4);
    e4s = computeE4s(n4e);
    epsilon = 1e-6;
    alpha = 1;
    gamma = 1;
    tau = 0.5;
    degreeF = 4;
    
    [STIMACR, MASSCR, intFSq4e, intFCR14e, intFCR24e, ...
     intFCR34e, gradsCR4e, s4e, dof, nrSides, area4e] ...
            = initPrimalDual(c4n, n4e, n4sDb, @f, degreeF);
    
    uCR = zeros(nrSides, 1);
    Lambda4e = zeros(size(n4e,1), 1);
%     uCR(dof) = ones(length(dof), 1);
%     
%     DuCR4e = computeGradientCR(uCR, gradsCR4e, s4e);
%     Lambda4e = DuCR4e./repmat(sqrt(sum(DuCR4e.^2,2)), 1, 2);
    
%     ind = (isnan(Lambda4e));
%     Lambda4e(ind) = 0;
    
    [uCR, ~] = solvePrimalDualFormulation(uCR, Lambda4e, alpha, epsilon, tau, ...
                                STIMACR, MASSCR, intFCR14e, intFCR24e, intFCR34e,...
                                gradsCR4e, s4e, dof, nrSides, area4e);
    
    figure;
    plotCR(c4n, n4e, uCR, 'Solution 1');
    
    [u,~,~,~] = ...
        tvRegPrimalDual(c4n,n4e,n4sDb,zeros(0,2),1/8,tau,3,epsilon,alpha,@f,uCR,Lambda4e);
    
    figure;
    plotCR(c4n, n4e, u, 'Solution 2');
end

function val = g(x,alpha,gamma)

    nP=size(x,1);
    val = zeros(nP,1);
    r=sqrt(sum(x.^2,2));

    temp = zeros(nP,1);
    temp(r<=1/6) = 1;
    r_temp = r(temp==1);
    val(temp==1) = alpha-12*(2-9*r_temp);

    temp = zeros(nP,1);
    temp(1/6<r & r<=1/3) = 1;
    r_temp = r(temp==1);
    val(temp==1) = alpha*(1+(6*r_temp-1).^gamma)-1./r_temp;

    temp = zeros(nP,1);
    temp(1/3<r & r<=1/2) = 1;
    r_temp = r(temp==1);
    val(temp==1) = 2*alpha+6*pi*sin(pi*(6*r_temp-2))-1./r_temp.*cos(pi*(6*r_temp-2));

    temp = zeros(nP,1);
    temp(1/2<r & r<=5/6) = 1;
    r_temp = r(temp==1);
    val(temp==1) = 2*alpha*(5/2-3*r_temp).^gamma+1./r_temp;

    temp = zeros(nP,1);
    temp(5/6<r & r<=1) = 1;
    r_temp = r(temp==1);
    val(temp==1) = -3*pi*sin(pi*(6*r_temp-5))+1./(2*r_temp).*(1+cos(pi*(6*r_temp-5)));

end

function value = f(x)
    value = ones(size(x,1),1);
end