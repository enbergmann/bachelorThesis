function energy = computeDiscreteEnergyCR(params, currData, v, vGradCR)
%% DOC
% Computes the discrete energy of the Crouzeix-Raviart function v
%   \alpha/2 ||v||^2_{L^2(Omega)} + ||\nabla_{CR} v||_{L^1(\Omega)}
%     - \int_\Omega fv dx
% for the nonconforming problem.
% 
% computeDiscreteEnergyCR.m
% input: params   - 'struct' with fields:
%                     parAlpha: 'double' containing the parameter alpha from
%                               the problem
%        currData - 'struct' with fields:
%                       area4e: areas for elements
%                     intRHS4s: '(nrSides x 1)-dimensional double array' where
%                               the j-th entry is the integral of f times the
%                               CR-basis function w.r.t. the j-th edge of the
%                               triangulation
%                       maMaCR: global CR mass matrix
%        v        - '(nrSides x 1)-dimensional double array' where the j-th row
%                   contains the coefficient of v w.r.t. the j-th edge of the
%                   triangulation 
%        vGradCR  - '(nrElems x 2)-dimensional double array' where the j-th row
%                   contains the gradient of v on the j-th triangle 
%
% output: energy - 'double' containing the discrete energy of u
  
%% INIT
  % extract necessary parameters from params
  parAlpha = params.parAlpha;

  % extract necessary information from currData
  area4e = currData.area4e;
  intRHS4s = currData.intRHS4s;
  maMaCR = currData.maMaCR;

%% MAIN
  energy = parAlpha/2*v'*maMaCR*v ...  % \alpha/2 ||v||^2_{L^2(Omega)}
    + area4e'*sqrt(sum(vGradCR.^2, 2)) ... % ||\nabla_{CR} v||_{L^1(\Omega)}
    - v'*intRHS4s; % \int_\Omega fv dx
end
