function y = f02(x, args)
  parAlpha = args(1, 1);
  parBeta = args(1, 2);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
 
  ind = 0<=r & r<=(1 - parBeta)/2;
  y(ind) = parAlpha - 12/(1 - parBeta)^2*r(ind) + 8/(1 - parBeta);

  ind = (1 - parBeta)/2<=r & r<=(1 + parBeta)/2;
  rTemp = r(ind);
  y(ind) = -parAlpha/parBeta*(rTemp - (1 + parBeta)/2) + 1./rTemp;

  ind = (1 + parBeta)/2<r & r<=1;
  rTemp = r(ind);
  y(ind) = -4/(parBeta - 1)^3*(16*rTemp.^2 - 9*(parBeta + 3)*rTemp + ...
    12*(parBeta + 1) - (3*parBeta + 1)./rTemp);
end
