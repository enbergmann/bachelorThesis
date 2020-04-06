function y = f02Gradient(x, args)
  parAlpha = args(1);
  
  r = sqrt(sum(x.^2, 2));

  y = zeros(length(r), 2);
  
  ind = 0<=r & r<=1/2;
  rTemp = r(ind);
  y(ind, :) = -(2*parAlpha-1./rTemp.^2).*x(ind, :)./rTemp;
  
  ind = 1/2<r & r<=1;
  rTemp = r(ind);
  y(ind, :) = 2*(64*rTemp-54+2./rTemp.^2).*x(ind, :)./rTemp;
end
