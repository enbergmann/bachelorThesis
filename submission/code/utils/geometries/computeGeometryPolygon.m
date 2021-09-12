function [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(j)
%% DOC
% Computes the data [c4n, n4e, n4sDb, n4sNb] for a mesh of the polygon
% approximating the circle with radius 1. The more red-refinements j are used
% the closer the polygon resembles the circle visually.
% 
% computeGeometryPolygon.m
% input: j - number of red-refinements to produce the mesh
%
% output: c4n   - coordinates for nodes
%         n4e   - nodes for elements
%         n4sDb - nodes for Dirichlet boundary sides
%         n4sNb - nodes for Neumann boundary sides

%% INIT
  c4n = [0 0; 0 1; -1 0; 0 -1; 1 0];
  n4e = [2 3 1;3 4 1;4 5 1;5 2 1];
  n4sDb = n4e(:, [1, 2]);
  n4sNb = [];
  
%% MAIN
  for k = 1:j
    [c4n, n4e, n4sDb, n4sNb] = refineUniformRed(c4n, n4e, n4sDb, n4sNb);
    temp = unique(n4sDb);
    c4n(temp, :) = ...
      c4n(temp,: )./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
  end
end
