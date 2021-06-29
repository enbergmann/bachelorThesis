function y = bubbleGradient(x)
  y = pi*[cos(pi*x(:, 1)).*sin(pi*x(:, 2)), sin(pi*x(:, 1)).*cos(pi*x(:, 2))];
end
