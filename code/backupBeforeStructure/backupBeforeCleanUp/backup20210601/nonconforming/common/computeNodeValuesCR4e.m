function nodeValues4e = computeNodeValuesCR4e(s4e, u)
%% DOC
% Computes the values of a CR function u in the three nodes of each triangle.
%
% computeNodeValuesCR4e.m
% input: s4e - sides for elements
%        u   - '(nrSides x 1)-dimensional double array' where the j-th row
%              contains the coefficient of the CR function u wrt. the j-th side
%              of the triangulation
%
% output: nodeValues4e - '(nrElems x 3)-dimensional double array' where the
%                        k-th entry of the j-th row contains the value of u in
%                        the k-th local node of the j-th triangle of the mesh

%% MAIN
  u = reshape(u, [length(u), 1]);
    % guarantee that u is a column
  nodeValues4e = [u(s4e(:, 3)) + u(s4e(:, 1)) - u(s4e(:, 2)), ... 
                  u(s4e(:, 1)) + u(s4e(:, 2)) - u(s4e(:, 3)), ...
                  u(s4e(:, 2)) + u(s4e(:, 3)) - u(s4e(:, 1))]; 
    % cf. thesis
end
