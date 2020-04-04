function y = f02ExactSolution(x)
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
  
  ind = r<=1/2;
  y(ind) = -2*r(ind) + 1;
end
