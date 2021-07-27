function plotAxisCR(c4n, n4e, u)
%% DOC
% Plots a CR function along the x- and y-axis on the mesh defined by [c4n,
% n4e].
%
% plotAxisCR.m
% input: c4n - coordinates for nodes
%        n4e - nodes for elements
%        u   - '(nrSides x 1)-dimensional double array' where the j-th row
%              contains the coefficient of the CR function w.r.t. the j-th side
%              of the triangulation

%% INIT
n4s = computeN4s(n4e);
nrSides = size(n4s, 1);
mid4s = computeMid4s(c4n, n4s);

temp1 = zeros(nrSides, 1);
temp2 = temp1;

%% MAIN
temp1(mid4s(:, 1)==0) = 1;
temp2(mid4s(:, 2)==0) = 1;

[y, I1] = sort(mid4s(temp1==1, 2));
[x, I2] = sort(mid4s(temp2==1, 1));

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
