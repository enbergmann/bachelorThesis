function plotGrayscale(c4n, n4e, x)
%% DOC
% Draws a grayscale picture on the triangulation defined by [c4n, n4e] with
% colors defined by the values in x. Call with 
%   x = mean(u(s4e), 2) 
% to draw a Crouzeix-Raviart function with coefficients u on the mesh [c4n,
% n4e] as a grayscale picture.
% 
% plotGrayscale.m
% input: c4n - coordinates for nodes
%        n4e - nodes for elements
%        x   - '(nrElems x 1)-dimensional double array' where the j-th entry
%              contains the value defining the color on the j-th triangle

%% INIT
  X1 = c4n(n4e(:, 1), 1);
  X2 = c4n(n4e(:, 2), 1);
  X3 = c4n(n4e(:, 3), 1);
  Y1 = c4n(n4e(:, 1), 2);
  Y2 = c4n(n4e(:, 2), 2);
  Y3 = c4n(n4e(:, 3), 2);
  X = [X1'; X2'; X3'];
  Y = [Y1'; Y2'; Y3'];
  
%% MAIN
  axis image;
  colormap gray;
  patch(X, Y, x', 'EdgeColor', 'none');
  drawnow;
end
