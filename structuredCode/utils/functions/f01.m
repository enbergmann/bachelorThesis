function y = f01(x, args)
  parAlpha = args(1, 1);
  parBeta = args(1, 2);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
 
  ind = r<=1/6;
  y(ind) = parAlpha - 12*(2-9*r(ind));
  
  ind = 1/6<r & r<=1/3;
  rTemp = r(ind);
  y(ind) = parAlpha*(1+(6*rTemp-1).^parBeta) - 1./rTemp;
  
  ind = 1/3<r & r<=1/2;
  rTemp = r(ind);
  arg = pi*(6*rTemp-2);
  y(ind) = 2*parAlpha + 6*pi*sin(arg) - 1./rTemp.*cos(arg);

  ind = 1/2<r & r<=5/6;
  rTemp = r(ind);
  y(ind) = 2*parAlpha*(5/2-3*rTemp).^parBeta + 1./rTemp;
  % TODO fix typo (+1/3 instead of +1/r) in thesis tex
  
  ind = 5/6<r & r<=1;
  rTemp = r(ind);
  arg = pi*(6*rTemp-5);
  y(ind) = -3*pi*sin(arg) + 1./(2*rTemp).*(1+cos(arg));
end
