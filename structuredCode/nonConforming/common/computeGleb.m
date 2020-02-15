function gleb = computeGleb(params, currData, ...
  output)
% TODO PastedGraphic-1.tiff 
% implement
% E_{NC}(u_{CR}) - \kappa_{CR}/\alpha ||h_\mathcal{T}(f-\alpha u_{CR})|| |f|_{1,2}
% where \kappa_{CR} is the constant from I_{NC} operator (GNUMO S 102, 
% \sqrt{1/48+1/j_{1,1}^2}\leq 0.298217419)
% TODO use gradient for f only when exact energy is known, hence make it 
% dependent on this to and then gradF is in params

  f = params.f;
  gradF = params.gradF;
  parAlpha = params.parAlpha;
  
  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  int1RHS4e = currData.int1RHS4e;
  int2RHS4e = currData.int2RHS4e;
  int3RHS4e = currData.int3RHS4e;
  length4s = currData.length4s;
  
  u = output.u;
  discreteEnergy = output.energyVec(end);
  
  kappaCR = 0.298217419;
  
  %TODO everything that was also computed in estimateError should prob be 
  % computed in maybe computeIntegralsWithF (rename it) and saved in currData
  % TODO ||f-alpha u||^2 was complely computed in estimateError, which must be
  % useful for ||f-alpha u|| here somehow
  termFSquared = integrate(@(n4p, Gpts4p, Gpts4ref)(f(Gpts4p).^2), ...
    c4n, n4e, 20, area4e);
    % termFSquared(elem) = ||f||^2_{L^2(T)}
  termMixed = u(s4e(:, 1)).*int1RHS4e + u(s4e(:, 2)).*int2RHS4e ... 
    + u(s4e(:, 3)).*int3RHS4e;
    % termMixed(elem) = (f,u)_{L^2(T)}
  termU = 1/3*area4e.*sum(u(s4e).^2, 2);
    % termU(elem) = ||u||_{L^2(T)}

  termGradFSquared = sum(sum(integrate(@(n4p, Gpts4p, Gpts4ref)(gradF(Gpts4p).^2), ...
    c4n, n4e, 20, area4e)));
    % termGradFSquared = ||\nabla f||^2_{L^2(\Omega)}

  integrals4e = max(...
    [length4s(s4e(:, 1)), length4s(s4e(:, 2)), length4s(s4e(:, 3))], [], 2) ...
    .*(termFSquared - 2*parAlpha*termMixed + parAlpha^2*termU);
  gleb = discreteEnergy ...
    - kappaCR\parAlpha*sqrt(sum(integrals4e)*termGradFSquared);
end
