function gleb = computeGleb(params, currData, ...
  output)
% TODO PastedGraphic-1.tiff 
% implement
% E_{NC}(u_{CR}) - \kappa_{CR}/\alpha ||h_\mathcal{T}(f-\alpha u_{CR})|| |f|_{1,2}
% where \kappa_{CR} is the constant from I_{NC} operator (GNUMO S 102, 
% \sqrt{1/48+1/j_{1,1}^2}\leq 0.298217419)
% TODO use gradient for f only when exact energy is known, hence make it 
% dependent on this to and then gradF is in params

  gradF = params.gradF;
  parAlpha = params.parAlpha;
  
  c4n = currData.c4n;
  n4e = currData.n4e;
  area4e = currData.area4e;
  s4e = currData.s4e;
  length4s = currData.length4s;
  
  discreteEnergy = output.energyVec(end);
  normOfDifference4e = output.normOfDifference4e;
  
  kappaCR = 0.298217419;
  
  termGradFSquared = sum(sum(integrate(@(n4p, Gpts4p, Gpts4ref)(gradF(Gpts4p).^2), ...
    c4n, n4e, 20, area4e)));
    % termGradFSquared = ||\nabla f||^2_{L^2(\Omega)}

  integrals4e = max(...
    [length4s(s4e(:, 1)), length4s(s4e(:, 2)), length4s(s4e(:, 3))], [], 2) ...
    .*normOfDifference4e;

  gleb = discreteEnergy ...
    - kappaCR\parAlpha*sqrt(sum(integrals4e)*termGradFSquared);
end
