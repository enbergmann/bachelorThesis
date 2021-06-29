function normDiffRhsSolCrSquared4e = ...
  computeNormDiffRhsSolCrSquared4e(params, currData, output)
%% DOC
% Computes 
%   ||f-\alpha u_{CR}||^2_{L^2(T)}
% for all T\in\mathcal{T} of a triangulation given by [c4n, n4e].
%
% computeNormOfDifference4e.m
% input: params   - 'struct' with fields:
%                                    f: 'function_handle' of the right-hand
%                                       side f
%                             parAlpha: 'double' containing the parameter alpha
%                                       from the problem
%                     degree4Integrate: 'uint64' up to which the integration in
%                                       integrate must be exact
%
%        currData - 'struct' with fields:
%                           c4n: coordinates for nodes
%                           n4e: nodes for elements
%                        area4e: areas for elements
%                           s4e: sides for elements
%                     int1RHS4e: '(nrElems x 1)-dimensional double array' where
%                                the j-th entry is the integral over the j-th
%                                element of f times the first local CR-basis
%                                function, i.e. the CR-basis function w.r.t.
%                                the first local edge of the j-th triangle
%                     int2RHS4e: as int1RHS4e for the second local CR-basis
%                                function
%                     int3RHS4e: as int1RHS4e for the third local CR-basis
%                                function
%
%        output   - 'struct' with fields:
%                      u: '(nrSides x 1)-dimensional double array' where the
%                         j-th row contains the coefficient of the CR solution
%                         u of the iteration w.r.t. the j-th side of the
%                         triangulation
%
% output: normDiffRhsSolCrSquared4e - '(nrElems x 1)-dimensional double array'
%                                     where the j-th row contains the integral
%                                     over the j-th triangle of the
%                                     triangulation of the right-hand side f
%                                     minus parAlpha times the CR solution u of
%                                     the iteration

%% INIT
  % extract necessary parameters from params
  f = params.f;
  parAlpha = params.parAlpha;
  degree4Integrate = params.degree4Integrate;

  % extract necessary information from currData
  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  int1RHS4e = currData.int1RHS4e;
  int2RHS4e = currData.int2RHS4e;
  int3RHS4e = currData.int3RHS4e;

  % extract necessary information from output
  u = output.u;

%% MAIN
  termRhsSquared4e = integrate(@(n4p, Gpts4p, Gpts4ref)(f(Gpts4p).^2), ...
    c4n, n4e, degree4Integrate, area4e);
    % termRhsSquare4e(elem) = ||f||^2_{L^2(T)}
  termMixed4e = u(s4e(:, 1)).*int1RHS4e + u(s4e(:, 2)).*int2RHS4e ... 
    + u(s4e(:, 3)).*int3RHS4e;
    % termMixed4e(elem) = (f, u)_{L^2(T)}
  termSolCrSquared4e = 1/3*area4e.*sum(u(s4e).^2, 2);
    % termSolCrSquared4e(elem) = ||u||^2_{L^2(T)}
    % NOTE u'*maMaCR*u would be faster to compute ||u|^2_{L^2(\Omega)} but 
    % the element-wise norm is necessary for computeGleb.m

  normDiffRhsSolCrSquared4e = termRhsSquared4e - 2*parAlpha*termMixed4e + ...
    parAlpha^2*termSolCrSquared4e;
end
