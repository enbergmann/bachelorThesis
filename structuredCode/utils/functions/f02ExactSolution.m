function y = f02ExactSolution(x, args)
  parBeta = args(1);

  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
  
  y(0<=r & r<=(1 - parBeta)/2) = 1;

  ind = (1 - parBeta)/2<=r & r<=(1 + parBeta)/2;
  y(ind) = -1/parBeta*r(ind) + (1 + parBeta)/(2*parBeta);

  y((1 + parBeta)/2<=r & r<=1) = 0;
end
