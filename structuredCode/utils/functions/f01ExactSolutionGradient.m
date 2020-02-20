function y = f01ExactSolutionGradient(x, args)
  parBeta = args(1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 2);
  
  ind = 1/6<r & r<=1/3;
  rTemp = r(ind);
  y(ind, :) = 6*parBeta*(6*rTemp-1).^(parBeta-1).*x(ind, :)./rTemp;
  
  ind = 1/2<=r & r<5/6;
  rTemp = r(ind);
  y(ind, :) = -6*parBeta*(5/2-3*rTemp).^(parBeta-1).*x(ind, :)./rTemp;
end

