function y = middleCircle(x, args)
  parAlpha = args(1, 1);
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);

  y(r<=1/2) = parAlpha;
end
