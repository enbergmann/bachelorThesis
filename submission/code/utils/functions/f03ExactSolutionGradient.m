function y = f03ExactSolutionGradient(x)
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 2);
  
  ind = 1/3<r & r<=2/3;
  rTemp = r(ind);
  y(ind, :) = (162*rTemp.^2 - 162*rTemp + 36).*x(ind, :)./rTemp;
end

