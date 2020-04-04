function y = f01Discontinuous(x, args)
  parAlpha = args(1, 1);
  parBeta = args(1, 2);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
 
  y(r<=1/6) = parAlpha;
  
  ind = 1/6<r & r<=1/3;
  rTemp = r(ind);
  y(ind) = parAlpha*(1+(6*rTemp-1).^parBeta) - 1./rTemp;
  
  y(1/3<r & r<=1/2) = 2*parAlpha;

  ind = 1/2<r & r<=5/6;
  rTemp = r(ind);
  y(ind) = 2*parAlpha*(5/2-3*rTemp).^parBeta + 1./rTemp;
  % TODO fix typo (+1/3 instead of +1/r) in thesis tex
end
