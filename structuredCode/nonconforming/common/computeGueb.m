function gueb = computeGueb(params, currData, output, outputLvlEnergy)
%% DOC TODO after questions are anwered
% Computes the guaranteed lower energy bound 
%   gleb = E_{NC}(u_{CR}) - 
%            \kappa_{CR}/\alpha ||h_\mathcal{T}(f-\alpha u_{CR})|| |f|_{1,2}
% w.r.t. the result of solvePrimalDualFormulation.m on the triangulation given 
% by [c4n, n4e] where
%   \kappa_{CR} = \sqrt{1/48 + 1/j_{1, 1}^2}\leq 0.298217419 
% (j_{1, 1} smallest positive root of the Bessel function of the first kind) is
% the constant in the L^2 approximation property of the I_{NC} operator (GNUMO
% p 102). 
%
% computeGueb.m
% input: params   - 'struct' with fields:
%                                gradF: 'function_handle' of the gradient of
%                                       the input signal f
%                             parAlpha: 'double' containing the parameter alpha
%                                       from the problem
%                     degree4Integrate: 'uint64' up to which the integration in
%                                       integrate must be exact
%
%        currData - 'struct' with fields:
%                          c4n: coordinates for nodes
%                          n4e: nodes for elements 
%                          s4e: sides for elements
%                       area4e: areas for elements
%                     length4s: lengths for sides
%
%        output   - 'struct' with fields: 
%                                     energyVec: '(nrIterations of current
%                                                level x 1)-dimensional double
%                                                array' where the j-th row
%                                                contains the discrete energy
%                                                of the j-th iterate of the
%                                                iteration on the current level
%                     normDiffRhsSolCrSquared4e: '(nrElems x 1)-dimensional
%                                                double array' where the j-th
%                                                row contains the integral over
%                                                the j-th triangle of the
%                                                triangulation of the square of
%                                                the difference of the
%                                                input signal f and the
%                                                product of parAlpha and the CR
%                                                solution u of the iteration
%
% output: gueb - 'double' containing the guaranteed upper energy bound

uCR = output.u;
maMaCR = currData.maMaCR;   
