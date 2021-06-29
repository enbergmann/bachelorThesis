function y = f02Gradient(x, args)
  parAlpha = args(1);
  parBeta = args(2);
  
  r = sqrt(sum(x.^2, 2));

  y = zeros(length(r), 2);
  
  ind = 0<r & r<=(1 - parBeta)/2;
  y(ind, :) = -12/(1 - parBeta)^2.*x(ind, :)./r(ind);
  
  ind = (1 - parBeta)/2<r & r<=(1 + parBeta)/2;
  rTemp = r(ind);
  y(ind, :) = (-parAlpha/parBeta - 1./rTemp.^2).*x(ind, :)./rTemp;

  ind = (1 + parBeta)/2<r & r<=1;
  rTemp = r(ind);
  y(ind, :) = -4/(parBeta - 1)^3*(32*rTemp - 9*(parBeta + 3) + ...
    (3*parBeta + 1)./rTemp.^2).*x(ind, :)./rTemp;
end
