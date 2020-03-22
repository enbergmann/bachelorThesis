function gleb = computeGleb(params, currData, output)
%% DOC
% Computes the guaranteed lower energy bound 
%   gleb = E_{NC}(u_{CR}) - 
%            \kappa_{CR}/\alpha ||h_\mathcal{T}(f-\alpha u_{CR})|| |f|_{1,2}
% w.r.t. the result of solvePrimalDualFormulation.m  on the triangulation given 
% by [c4n, n4e] where
%   \kappa_{CR} = \sqrt{1/48 + 1/j_{1, 1}^2}\leq 0.298217419 
% (j_{1, 1} smallest positive root of the Bessel function of the first kind) is
% the constant in the L^2 approximation property of the I_{NC} operator (GNUMO
% p 102). 
%
% computeGleb.m
% input:  params   - 'struct' with fields:
%                                 gradF: 'function_handle' of the gradient of
%                                        the right-hand side f
%                              parAlpha: 'double' containing the parameter 
%                                        alpha from the problem
%                      degree4Integrate: 'uint64' up to which the integration
%                                        in integrate must be exact
%
%         currData - 'struct' with fields:
%                           c4n: coordinates for nodes
%                           n4e: nodes for elements 
%                           s4e: sides for elements
%                        area4e: areas for elements
%                      length4s: lengths for sides
%
%         output   - 'struct' with fields: 
%                               energyVec: '(nrIterations of current level 
%                                          x 1)-dimensional double array' where
%                                          the j-th row contains the discrete 
%                                          energy of the j-th iterate of the 
%                                          iteration on the current level
%                      normOfDifference4e: TODO
%
% output: gleb     - 'double' containing the guaranteed lower energy bound

%% INIT
  % extract necessary parameters from params
  gradF = params.gradF;
  parAlpha = params.parAlpha;
  degree4Integrate = params.degree4Integrate;
  
  % extract necessary information from currData
  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  length4s = currData.length4s;
  
  % extract necessary information from output
  discreteEnergy = output.energyVec(end);
  normOfDifference4e = output.normOfDifference4e;
  
  % define remaining parameters
  kappaCR = 0.298217419;
  
%% MAIN
  termGradFSquared = sum(...
    integrate(@(n4p, Gpts4p, Gpts4ref)(sum(gradF(Gpts4p).^2, 2)), ...
    c4n, n4e, degree4Integrate, area4e));
    % ||\nabla f||^2_{L^2(\Omega)}

  termScaledDifference = sum(...
    max([length4s(s4e(:, 1)), length4s(s4e(:, 2)), length4s(s4e(:, 3))], ...
    [], 2).*normOfDifference4e);
    % ||h_\\mathcal{T}(f - \alpha u_{CR})||^2_{L^2(\Omega)}

  % compute gleb
  gleb = discreteEnergy - ...
    kappaCR/parAlpha*sqrt(termScaledDifference*termGradFSquared);
end
