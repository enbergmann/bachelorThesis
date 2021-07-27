function y = f04Gradient(x, args)
  parAlpha = args(1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 2);
  
  ind = 0<r & r<=1/3;
  rTemp = r(ind);
  y(ind, :) = (34992*rTemp.^3 - 18225*rTemp.^2 + 2160*rTemp).*x(ind, :)./rTemp;
  
  ind = 1/3<r & r<=2/3;
  rTemp = r(ind);
  y(ind, :) = (parAlpha*(162*rTemp.^2 - 162*rTemp + 36) - ...
    1./rTemp.^2).*x(ind, :)./rTemp;

  ind = 2/3<r & r<=1;
  rTemp = r(ind);
  y(ind, :) = (3645*rTemp.^2 - 6048*rTemp + 2592 - ...
    81./rTemp.^2).*x(ind, :)./rTemp;
end
