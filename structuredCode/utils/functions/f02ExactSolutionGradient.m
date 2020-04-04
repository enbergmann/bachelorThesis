function y = f02ExactSolutionGradient(x)
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 2);
  
  ind = r<=1/2;
  y(ind, :) = -2*x(ind, :)./r(ind);
end

