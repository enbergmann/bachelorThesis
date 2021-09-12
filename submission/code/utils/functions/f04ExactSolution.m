function y = f04ExactSolution(x)
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
  
  y(r<=1/3) = 1;

  ind = 1/3<r & r<=2/3;
  rTemp = r(ind);
  y(ind) = 54*rTemp.^3 - 81*rTemp.^2 + 36*rTemp - 4;

  y(2/3<r & r<=1) = 0;
end
