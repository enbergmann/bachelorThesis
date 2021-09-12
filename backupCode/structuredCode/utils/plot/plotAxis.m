function plotAxis(c4n, u)
%% DOC
% Plots a function along the x- and y-axis on the mesh defined by c4n.
%
% plotAxis.m
% input: c4n - coordinates for nodes of the coarse mesh
%        u   - '(nrNodes x 1)-dimensional double array' where the j-th row
%              contains the coefficient of the function w.r.t. the j-th node
%              of the triangulation

%% INIT
temp1 = zeros(size(c4n, 1), 1);
temp2 = temp1;

temp1(c4n(:, 1)==0) = 1;
temp2(c4n(:, 2)==0) = 1;

%% MAIN
[y, I1] = sort(c4n(temp1==1, 2));
[x, I2] = sort(c4n(temp2==1, 1));

u1 = u(temp1==1);
u2 = u(temp2==1);

u1 = u1(I1);
u2 = u2(I2);

if length(u1)<1e2
  plot(y, u1, '-*r', 'LineWidth', 3);
  hold on
  plot(x, u2, '-ob', 'LineWidth', 1);
else
  plot(y, u1, '-r', 'LineWidth', 3);
  hold on
  plot(x, u2, '-b', 'LineWidth', 1);
end
