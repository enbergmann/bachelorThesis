function E = energyAnalyticExample(red,alpha,delta)
  u = @(x) gUexact(x,alpha,delta);
  [c4n,n4e] = computeGeometryPolygon(red);
  uNC = interpolationNC(u,c4n,n4e,computeN4s(n4e));
  E = exactDiscreteEnergy(c4n,n4e,computeS4e(n4e),uNC);
end
