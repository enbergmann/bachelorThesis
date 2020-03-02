function normOfDifference4e = ...
  computeNormOfDifference4e(params, currData, output)

  f = params.f;
  parAlpha = params.parAlpha;

  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  int1RHS4e = currData.int1RHS4e;
  int2RHS4e = currData.int2RHS4e;
  int3RHS4e = currData.int3RHS4e;
    % intjRHS4e(elem) = \int_T psi_j *f dx

  u = output.u;

  termFSquared = integrate(@(n4p, Gpts4p, Gpts4ref)(f(Gpts4p).^2), ...
    c4n, n4e, 20, area4e);
    % termFSquared(elem) = ||f||^2_{L^2(T)}
  termMixed = u(s4e(:, 1)).*int1RHS4e + u(s4e(:, 2)).*int2RHS4e ... 
    + u(s4e(:, 3)).*int3RHS4e;
    % termMixed(elem) = (f,u)_{L^2(T)}
  termU = 1/3*area4e.*sum(u(s4e).^2, 2);
    % termU(elem) = ||u||_{L^2(T)}

  normOfDifference4e = (termFSquared - 2*parAlpha*termMixed + ...
    parAlpha^2*termU);
end
