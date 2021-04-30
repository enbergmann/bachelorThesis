function y = f03(x, args)
  parAlpha = args(1, 1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
 
  ind = r<=1/3;
  rTemp = r(ind);
  y(ind) = parAlpha + 8748*rTemp.^4 - 6075*rTemp.^3 + 1080*rTemp.^2;
  
  ind = 1/3<r & r<=2/3;
  rTemp = r(ind);
  y(ind) = parAlpha*(54*rTemp.^3 - 81*rTemp.^2 + 36*rTemp - 4) + 1./rTemp;
  
  ind = 2/3<r & r<=1;
  rTemp = r(ind);
  y(ind) = 1215*rTemp.^3 - 3024*rTemp.^2 + 2592*rTemp - 864 + 81./rTemp;
end
