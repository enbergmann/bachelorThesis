function y = f03(x, args)
  parAlpha = args(1, 1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
 
  ind = r<=1/3;
  rTemp = r(ind);
  y(ind) = parAlpha - 216*rTemp.^2 + 81*rTemp;
  
  ind = 1/3<r & r<=2/3;
  rTemp = r(ind);
  y(ind) = parAlpha*(54*rTemp.^3 - 81*rTemp.^2 + 36*rTemp - 4) + 1./rTemp;
  
  ind = 2/3<r & r<=1;
  rTemp = r(ind);
  y(ind) = 216*rTemp.^2 - 405*rTemp + 216 - 27./rTemp;
end
