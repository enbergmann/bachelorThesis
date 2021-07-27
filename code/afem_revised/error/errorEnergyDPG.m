function  err  = errorEnergyDPG( gradF, u4e, c4n, n4e, kU, baseU, gradU, degreeF)
%ERRORENERGYDPG Computes the energy error err = |f-u| for the solution u to the 
%   primal dPG method and the exact solution f.
%   
% input:    f           - exact solution 
%           u4e         - solution to the Primal DPG method for the Poisson
%                         model problem for every triangle                      
%           c4n         - coordinates for nodes
%           n4e         - nodes for element
%           kU          - total polynomial degree of P_{k_u}(\mathcal{T})
%           baseU       - basis for P_{k_u}(\mathcal{T}) as a part of the
%                         trial space, for simplicity only on the reference
%                         triangle; order: base functions in vertices,
%                         on first side, second side, third side, inner
%                         points
%           nCPu        - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIPu        - number of base functions in inner points
%           degreeF     - polynomial degree for Gaussian quadrature rule
%                         for f
%
% output:   err         - the error |f-u| in the energy norm
%

% Copyright 2017 Enrico Bergmann, Leonard Richter-Matthies, Maximilian Schade
% reference version: MATLAB R2016a

%% INITIALIZATION
area4e = computeArea4e(c4n, n4e);
nElem = size(n4e, 1);
nBaseU = size(baseU, 2);
degree = max(kU,degreeF); % degree for gaussian integration

[~, gradU4e] = transformAffine(c4n, n4e, baseU, gradU);
    
%% COMPUTE ERROR
% integrate (gradF - gradU)^2 on each element
INT = @(n4p, Gpts4p, Gpts4ref)...
   sum((gradF(Gpts4p)-evalG(Gpts4ref, nElem, nBaseU, gradU4e, u4e)).^2, 2);

err = sum(integrate(INT, c4n, n4e, (degree-1)^2, area4e));
err = sqrt(err);
end

function sumsum = evalG(Gpts, nElem, nBaseU, grad, u4e)
% computes gradU on each element 
% (sums the base functions with the weight given by u4e)
sumsum = zeros(nElem,2);
for k=1:nBaseU
sumsum = sumsum + [u4e(:,k), u4e(:,k)].*grad{k}(Gpts);
end
end