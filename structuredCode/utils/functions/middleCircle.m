function y = middleCircle(x, args)
  parAlpha = args(1, 1);
  r = sqrt(sum(x.^2, 2));
  y = zeros(length(r), 1);

  % TODO
  y(-1/2<x(:, 1) & x(:, 1) < 1/2 & -1/2< x(:, 2) & x(:, 2)<1/2) = parAlpha;
end
