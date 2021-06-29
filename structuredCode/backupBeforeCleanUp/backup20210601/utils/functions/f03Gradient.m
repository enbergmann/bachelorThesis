function y = f03Gradient(x, args)
  parAlpha = args(1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 2);
  
  ind = 0<r & r<=1/3;
  rTemp = r(ind);
  y(ind, :) = (-432*rTemp + 81).*x(ind, :)./rTemp;
  
  ind = 1/3<r & r<=2/3;
  rTemp = r(ind);
  y(ind, :) = (parAlpha*(162*rTemp.^2 - 162*rTemp + 36) - ...
    1./rTemp.^2).*x(ind, :)./rTemp;

  ind = 2/3<r & r<=1;
  rTemp = r(ind);
  y(ind, :) = (432*rTemp - 405 + 27./rTemp.^2).*x(ind, :)./rTemp;
end
