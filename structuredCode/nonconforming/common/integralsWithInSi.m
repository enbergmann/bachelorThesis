function [int1InSi4e, int2InSi4e, int3InSi4e, intInSi4s] = ...
    integralsWithInSi(params, currData)
%% DOC
% Computes the element-wise integrals of a function f times the three local
% Crouzeix-Raviart basis functions on an triangle and the integrals of f times
% the global Crouzeix-Raviart basis functions, all with respect to the
% triangulation given by [c4n, n4e].
% 
% integralsWithInSi.m
% input: params   - 'struct' with fields:
%                                    f: 'function_handle' of the input signal f
%                     degree4Integrate: 'uint64' up to which the integration in
%                                       integrate must be exact
%        currData - 'struct' with fields:
%                         c4n: coordinates for nodes
%                         n4e: nodes for elements
%                         s4e: sides for elements
%                      area4e: areas for elements
%                     nrSides: number of sides
%                     nrElems: number of elements
%
% output: int1InSi4e - '(nrElems x 1)-dimensional double array' where the j-th
%                      entry is the integral over the j-th element of f times
%                      the first local CR-basis function, i.e. the CR-basis
%                      function w.r.t. the first local edge of the j-th
%                      triangle
%         int2InSi4e - as int1InSi4e for the second local CR-basis function
%         int3InSi4e - as int1InSi4e for the third local CR-basis function
%         intInSi4s  - '(nrSides x 1)-dimensional double array' where 
%                      the j-th entry is the integral of f times the CR-basis
%                      function w.r.t. the j-th edge

%% INIT
  % extract necessary parameters from params
  f = params.f;
  degree4Integrate = params.degree4Integrate;
  
  % extract necessary parameters from params
  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  nrSides = currData.nrSides;
  nrElems = currData.nrElems;
  s4e = currData.s4e;
  
  % compute barycentric coordinates
  % NOTE could be a separate function if necessary
  a1 = c4n(n4e(:, 1), :); 
  a2 = c4n(n4e(:, 2), :);
  a3 = c4n(n4e(:, 3), :);
    % the j-th row of ak contains the coordinates Pk of the k-th local node of
    % the j-th triangle T = conv({P1, P2, P3}), i.e. ak(j, :) = [Pk1 Pk2]

  matB4e = [a1 - a3, a2 - a3];
    % define matrix matB for affine transformation 
    %   B(x) = matB(x) + P3 
    % from the reference triangle Tref = conv({[1; 0], [0; 1], [0; 0]}) to a
    % triangle T = conv({P1, P2, P3}) where matB = [P1 - P3, P2 - P3] in row
    % notation, i.e. a matrix [a c; b d] is written as [a b c d]
    % this means B([1; 0]) = P1, B([0; 1]) = P2, B([0; 0]) = P3

  matBinv4e = ...
    1./(matB4e(:, 1).*matB4e(:, 4) - matB4e(:, 2).*matB4e(:, 3))...
    .*[matB4e(:, 4), -matB4e(:, 2), -matB4e(:, 3), matB4e(:, 1)];
    % inverse matrix matBinv of the affine transformation B utilizing the
    % formula for the inverse of a 2x2 matrix
    %   [a c; b d]^{-1} = 1/(ad - bc)*[d -c; -b a]
    % in row notation, i.e. the matrix [d -c; -b a] is written as [d -b -c a]
    % used to define barycentric coordinates by means of affine transformation
    %   Binv(y) = matBinv(y) - matBinv(P3) = matBinv(y - P3)
    % from a triangle T = conv({P1, P2, P3}) to the 
    % reference triangle Tref = conv({[1; 0], [0; 1], [0; 0]})
    % this means Binv(P1) = [1; 0], Binv(P1) = [0; 1], Binv(P3) = [0; 0],
  
  funcLambda1 = @(y)(sum([matBinv4e(:, 1) matBinv4e(:, 3)].*(y - a3), 2));
    % first component of Binv(y), i.e. funcLambda1(P1) = 1 and 
    % funcLambda1(P2) = funcLambda1(P3) = 0, hence, by affine linearity of
    % funcLambda1 it is the barycentric coordinate w.r.t. P1
  funcLambda2 = @(y)(sum([matBinv4e(:, 2) matBinv4e(:, 4)].*(y - a3), 2));
    % second component of Binv(y), i.e. funcLambda2(P2) = 1 and 
    % funcLambda2(P1) = funcLambda2(P3) = 0, hence, by affine linearity of
    % funcLambda2 it is the barycentric coordinate w.r.t. P2
  funcLambda3 = @(y)(1 - funcLambda1(y) - funcLambda2(y));
    % affine function with funcLambda3(P3) = 1 and 
    % funcLambda3(P1) = funcLambda3(P2) = 0, hence funcLambda3 is the 
    % barycentric coordinate w.r.t. P3

%% MAIN
  % compute integrals \int_T f*\psiLocal_j dx
  int1InSi4e = integrate(@(n4p, Gpts4p, Gpts4ref)(...
    f(Gpts4p).*(1 - 2*funcLambda3(Gpts4p))), ...
    c4n, n4e, degree4Integrate, area4e);
  int2InSi4e = integrate(@(n4p, Gpts4p, Gpts4ref)(...
    f(Gpts4p).*(1 - 2*funcLambda1(Gpts4p))), ...
    c4n, n4e, degree4Integrate, area4e);
  int3InSi4e = integrate(@(n4p, Gpts4p, Gpts4ref)(...
    f(Gpts4p).*(1 - 2*funcLambda2(Gpts4p))), ...
    c4n, n4e, degree4Integrate, area4e);

  % compute integrals \int_\Omega f*\psi_j dx
  intInSi4s = zeros(nrSides, 1); 
  for elem = 1 : nrElems
    intInSi4s(s4e(elem, :)) = intInSi4s(s4e(elem,:)) + ...
      [int1InSi4e(elem), int2InSi4e(elem), int3InSi4e(elem)]';
  end
end
