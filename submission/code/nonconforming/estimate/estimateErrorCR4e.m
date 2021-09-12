function [eta4e, etaVol4e, etaJumps4e] = ...
    estimateErrorCR4e(params, currData, output)
%% DOC
% Computes the refinement indicator \eta(T) = \etaVol(T) + \etaJumps(T) with
%   eta_{Vol}(T) = |T|^{2/d}||f - \alpha u_{CR}||^2_{L^2(T)} and
%   eta_{Jumps}(T) = ...
%     |T|^{\gamma/d}\sum_{F\in\mathal{F}(T)}||[u_{CR}]_F||_{L^1(F)}
% for the solution u_{CR} of the primal dual iteration and for all
% T\in\mathcal{T} w.r.t. the triangulation given by [c4n, n4e].
%
% estimateErrorCR4e.m
% input: params   - 'struct' with fields:
%                     parGamma: 'double' containing the parameter \gamma for
%                               the refinement indicator
%                            d: 'uint64' containing the dimension (only tested
%                                for d = 2)
%                                   
%        currData - 'struct' with fields:
%                             area4e: areas for elements
%                                s4e: sides for elements
%                     l1NormOfJump4s: '(nrSides x 1)-dimensional double array'
%                                     where the j-th row contains the L1 norm
%                                     of the jump of u across the j-th side of
%                                     the triangulation (if a side is an outer
%                                     edge the jump is defined as the
%                                     restriction of u to this side)
%
%        output   - 'struct' with fields: 
%                     normDiffInSiSolCrSquared4e: '(nrElems x 1)-dimensional
%                                                 double array' where the j-th
%                                                 row contains the integral
%                                                 over the j-th triangle of the
%                                                 triangulation of the square
%                                                 of the difference of the
%                                                 input signal f and the
%                                                 product of parAlpha and the
%                                                 CR solution u of the
%                                                 iteration
%
% output: eta4e      - '(nrElems x 1)-dimensional double array' where the j-th
%                      row contains the contribution from the refinement
%                      indicator \eta on the j-th triangle of the triangulation
%         etaVol4e   - '(nrElems x 1)-dimensional double array' where the j-th
%                      row contains the contribution from the volume part
%                      \eta_{Vol} of the refinement indicator on the j-th
%                      triangle of the triangulation
%         etaJumps4e - '(nrElems x 1)-dimensional double array' where the j-th
%                      row contains the contribution from the jump part
%                      \eta_{Jumps} of the refinement indicator on the j-th
%                      triangle of the triangulation

%% INIT
  % extract necessary parameters from params
  parGamma = params.parGamma;
  d = params.d;

  % extract necessary information from currData
  area4e = currData.area4e;
  s4e = currData.s4e;
  l1NormOfJump4s = currData.l1NormOfJump4s;

  % extract necessary information from output
  normDiffInSiSolCrSquared4e = output.normDiffInSiSolCrSquared4e;
  
%% MAIN
  etaVol4e = area4e.^(2/d).*normDiffInSiSolCrSquared4e;
  etaJumps4e = area4e.^(parGamma/d).*sum(l1NormOfJump4s(s4e), 2);

  eta4e = etaVol4e + etaJumps4e;
end
