function [eta4e, etaVol4e, etaJumps4e] = ...
    estimateErrorCR4e(params, currData, output)
%% DOC
% Computes the refinement indicator \eta(T) = \etaVol(T) + \etaJumps(T) with
%   eta_{Vol}(T) = |T|^{2/n}||f-\alpha u_{CR}||^2_{L^2(T)} and
%   eta_{Jumps}(T) = ...
%     |T|^{\beta/n}\sum_{F\in\mathal{F}(T)}||[u_{CR}]_F||_{L^1(F)}
% for all T\in\mathcal{T} w.r.t. the triangulation given by [c4n, n4e].
%
% estimateErrorCR4e.m
% input: params   - 'struct' with fields:
%                     beta4Estimate: 'double' containing the parameter \beta
%                                    from the problem
%                        n4Estimate: 'uint64' containing the dimension
%                                   
%        currData - 'struct' with fields:
%                          n4e: nodes for elements 
%                       area4e: areas for elements
%                          s4e: sides for elements
%                     length4s: lengths for sides
%                          e4s: elements for sides
%
%        output   - 'struct' with fields: 
%                                             u: '(nrSides x 1)-dimensional
%                                                double array' where the j-th
%                                                row contains the coefficient
%                                                of the solution u of the
%                                                iteration w.r.t. the j-th
%                                                edge of the triangulation 
%                     normDiffRhsSolCrSquared4e: '(nrElems x 1)-dimensional
%                                                double array' where the j-th
%                                                row contains the integral over
%                                                the j-th triangle of the
%                                                triangulation of the square of
%                                                the difference of the
%                                                right-hand side f and the
%                                                product of parAlpha and the CR
%                                                solution u of the iteration
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
  beta4Estimate = params.beta4Estimate;
  n4Estimate = params.n4Estimate;

  % extract necessary information from currData
  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  length4s = currData.length4s;
  e4s = currData.e4s;

  % extract necessary information from output
  u = output.u;
  normDiffRhsSolCrSquared4e = output.normDiffRhsSolCrSquared4e;
  
%% MAIN
  % compute \eta_{Vol}
  etaVol4e = area4e.^(2/n4Estimate).*normDiffRhsSolCrSquared4e;

  % compute \eta_{Jumps}
  nodeValues4e = computeNodeValuesCR4e(s4e, u); 
  absNodeJumps4s = computeAbsNodeJumps4s(n4e, e4s, nodeValues4e);

  termJumps4e = 1/4*(...
                length4s(s4e(:, 1)).*sum(absNodeJumps4s(s4e(:, 1), :), 2)...
                + length4s(s4e(:, 2)).*sum(absNodeJumps4s(s4e(:, 2), :), 2)...
                + length4s(s4e(:, 3)).*sum(absNodeJumps4s(s4e(:, 3), :), 2)); 
    % for any edge F := conv({a, b}) = Tp \cap Tm, m := (a+b)/2,
    % u := u_{CR}, up := u|_Tp, and um := u|_Tm, and d := |up - um| it holds
    %   ||[u]_F||_{L^1(F)} = \int_F |up - um|(x) dx = \int_a^b d(x) dx 
    %                      = \int_a^m d(x) dx + \int_m^b d(x) dx
    % by continuity of the CR function u in m and affine linearity of up and
    % um the difference d = |up - um| satisfies d(m) = 0 and is affine linear
    % on both [a, m] and [m, b] with
    %   d(x) = (d(m) - d(a))/(m - a) (x - a) + d(a) 
    %        = -d(a)/(m - a) (x - a) + d(a)             
    % on [a, m] and
    %   d(x) = (d(b) - d(m))/(b - m) (x - m) + d(m)
    %        = d(b)/(b - m) (x - m)
    % on [m, b]
    % the midpoints rule, which is exact for affine linear functions, implies
    %   ||[u]_F||_{L^1(F)} = |m - a|d((m + a)/2) + |b - m|d((b + m)/2)
    % by definition of m it holds m - a = (b - a)/2 = b - m, hence
    %   ||[u]_F||_{L^1(F)} = |b - a|/2(d((m + a)/2) + d((b + m)/2))
    %                      = |b - a|/2(-d(a)/(m - a) ((m + a)/2 - a) + d(a) 
    %                                  + d(b)/(b - m) ((b + m)/2 - m))
    %                      = |b - a|/2(-d(a)/(m - a) (m - a)/2 + d(a) 
    %                                  + d(b)/(b - m) (b - m)/2)
    %                      = |b - a|/2(-d(a)/2 + d(a) + d(b)/2)
    %                      = |b - a|/2(d(a)/2 + d(b)/2)
    %                      = |b - a|/4(d(a) + d(b))
    % altogether
    %   ||[u]_F||_{L^1(F)} = |b - a|/4(|up(b) - um(b)| + |up(a) - um(a)|)

  etaJumps4e = area4e.^(beta4Estimate/n4Estimate).*termJumps4e;

  % compute \eta
  eta4e = etaVol4e + etaJumps4e;
end
