function y = f01ExactSolution(x, args)
  parBeta = args(1);
  
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);
  
  y(r<=1/6) = 1;
  
  ind = 1/6<r & r<=1/3;
  y(ind) = 1 + (6*r(ind) - 1).^parBeta;
  
  y(1/3<r & r<=1/2) = 2;
  
  ind = 1/2<r & r<=5/6;
  y(ind) = 2*(5/2 - 3*r(ind)).^parBeta;
end
