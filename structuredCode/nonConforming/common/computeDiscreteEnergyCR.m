% TODO recalculate this stuff (TIens looks different somehow)
%
function E = computeDiscreteEnergyCR(params, currData, u, gradCRu)
% Compute the discrete energy for the nonconforming problem of the
% Crouzeix-Raviart function u. 
% 
% computeDiscreteEnergyCR.m
% input:  params   - 'struct' with fields:
%                       parAlpha: parameter alpha from the problem
%         currData - 'struct' with fields:
%                         area4e: areas for elements
%                       intRHS4s: integral of f times the CR-basis function
%                                 wrt. j-th side in the j-th component
%                         maMaCR: global CR mass matrix
%
% output: E        - 'double' containing the discrete energy of u
  
  % extract necessary data
  parAlpha = params.parAlpha;

  area4e = currData.area4e;
  intRHS4s = currData.intRHS4s;
  maMaCR = currData.maMaCR;

  % compute energy
  E = area4e'*sqrt(sum(gradCRu.^2, 2)) ... % L^1 norm uf gradCR u
               - u'*intRHS4s ... % integral f u
               + parAlpha/2*u'*maMaCR*u;    % L^2 norm uf u
end
