function y = f01(x, args)
  parAlpha = args(1, 1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
 
  ind = r<=1/2;
  rTemp = r(ind);
  y(ind) = parAlpha*(-2*rTemp + 1) + 1./rTemp;

  ind = 1/2<r & r<=1;
  rTemp = r(ind);
  y(ind) = 64*rTemp.^2 - 108*rTemp + 48 - 4./rTemp;
end
