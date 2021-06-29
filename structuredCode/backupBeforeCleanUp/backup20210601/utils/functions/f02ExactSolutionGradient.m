function y = f02ExactSolutionGradient(x, args)
  parBeta = args(1);

  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 2);
  
  ind = (1 - parBeta)/2<r & r<=(1 + parBeta)/2;
  y(ind, :) = -1/parBeta*x(ind, :)./r(ind);
end

