function [int1RHS4e, int2RHS4e, int3RHS4e, intRHS4s] = ...
    integralsWithF4e(currData, f, degree)
% Computes the integrals of a function f times a Crouzeix-Raviart-basis
% function, both global and local, with respect to the triangulation given by
% [c4n, n4e].
% 
% integralsWithF4e.m
% input:  currData   - 'struct' with fields:
%                           c4n: coordinates for nodes
%                           n4e: nodes for elements
%                        area4e: area for elements
%         degree     - 'uint64' up to which the integration in integrate must 
%                      be exact
%         f          - 'function_handle' of the function f in the integral
%
% output: int1RHS4e  - '(nrElems x 1)-dimensional double array' where 
%                      the j-th entry is the integral over the j-th element of
%                      f times the first local CR-basis function, i.e. the
%                      CR-basis function wrt. the first local edge of the j-th
%                      triangle
%         int2RHS4e  - as int1RHS4e for the second local CR-basis function
%         int3RHS4e  - as int1RHS4e for the third local CR-basis function
%         intRHS4s   - '(nrSides x 1)-dimensional double array' where 
%                      the j-th entry is the integral of f times the CR-basis
%                      function wrt. j-th edge

  % extract necessary dCR-basis function wrt. the first local edge of the j-thata
  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  nrSides = currData.nrSides;
  nrElems = currData.nrElems;
  s4e = currData.s4e;

  %TODO might just replace that with defining barCoords on Tref and using
  %Gpts4ref for integrate 
  
  % compute barycentric coordinates

  % vertices
  a1 = c4n(n4e(:, 1), :);
  a2 = c4n(n4e(:, 2), :);
  a3 = c4n(n4e(:, 3), :);

  % define barycentric coordinates by means of affine transformation
  % varPhi(x)+a3 where varPhi = [a1-a3,a2-a3] '2x2-matrix' in row notation
  % where [a c; b d] --> [a b c d]
  varPhi = [a1-a3, a2-a3];
  % compute Psi=Phi_inverse
  varPsi = repmat(1./(varPhi(:,1).*varPhi(:,4) - varPhi(:,2).*varPhi(:,3)),...
    1, 4).*[varPhi(:,4), -varPhi(:,2), -varPhi(:,3), varPhi(:,1)];
  % define barycentric coordinates
  varLambda1 = @(z)(varPsi(:, 1).*(z(:, 1)-a3(:, 1)) + varPsi(:, 3).*...
                                           (z(:, 2)-a3(:, 2)));
  varLambda2 = @(z)(varPsi(:, 2).*(z(:, 1)-a3(:, 1)) + varPsi(:, 4).*...
                                           (z(:, 2)-a3(:, 2)));
  varLambda3 = @(z)(1 - varLambda1(z) - varLambda2(z));

  % compute integrals \int_\Omega f*\psi_j dx
  int1RHS4e = integrate(@(n4p, Gpts4p, Gpts4ref)(...
    f(Gpts4p).*(1 - 2*varLambda3(Gpts4p))), c4n, n4e, degree+1, area4e);
  int2RHS4e = integrate(@(n4p, Gpts4p, Gpts4ref)(...
    f(Gpts4p).*(1 - 2*varLambda1(Gpts4p))), c4n, n4e, degree+1, area4e);
  int3RHS4e = integrate(@(n4p, Gpts4p, Gpts4ref)(...
    f(Gpts4p).*(1 - 2*varLambda2(Gpts4p))), c4n, n4e, degree+1, area4e);
  
  % compute integrals
  intRHS4s = zeros(nrSides, 1); % int_\Omega f*u dx
  for elem = 1 : nrElems
    intRHS4s(s4e(elem, :)) = intRHS4s(s4e(elem,:)) + ...
      [int1RHS4e(elem), int2RHS4e(elem), int3RHS4e(elem)]';
  end
end
