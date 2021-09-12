function y = middleSquare(x)
  y = zeros(size(x, 1), 1);
  y(-1/2<x(:, 1) & x(:, 1) < 1/2 & -1/2< x(:, 2) & x(:, 2)<1/2) = 100;
end
