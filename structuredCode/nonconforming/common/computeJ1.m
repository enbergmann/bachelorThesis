function vJ1 = prolongationJ1(c4n, n4e, n4sB, v)
%% DOC
% Computes the enriching operator of a CR function v defined by 
%   J_1 v (z) := |\Tcal(z)|^{-1} \sum_{T\in\Tcal(z)}v|_T(z) 
% for all z\in\Ncal(\Omega) as a P_1(\Tcal)\cap C_0(\overline\Omega) (Courant)
% function with respect to the mesh defined by [c4n, n4e].
%
% prolongationJ1.m
% input: c4n    - coordinates for nodes of the coarse mesh
%        n4e    - nodes for elements of the coarse mesh
%        n4sB   - nodes for boundary sides of the coarse mesh
%        v      - '(nrSides x 1)-dimensional double array' where the j-th row
%                   contains the coefficient of the CR function w.r.t. the j-th
%                   side of the triangulation
%
% output: vJ1   - '(nrNodes x 1)-dimensional double array' where the j-th row
%                  contains the coefficient of J1v w.r.t. the j-th node of the
%                  triangulation

%% INIT
  s4e = computeS4e(n4e); 

  nodeValuesCR4e = computeNodeValuesCR4e(s4e, v);
  innerNodes = setdiff(1:size(c4n, 1), unique(n4sB(:)));
  % TODO test this below, should be equivalent and hence c4n wouldnt be needed
  % innerNodes = setdiff(1:max(max(n4e)), unique(n4sB(:)));

%% MAIN
  % compute enriching operator
  nodeValuesJ1 = zeros(size(nodeValuesCR4e));
  for j = innerNodes
    indNode = find(n4e == j);
    localVals4node = nodeValuesCR4e(indNode);
    nodeValuesJ1(indNode) = sum(localVals4node)/length(localVals4node);
  end
  % vJ1 = ...;
  %  TODO not done here, need to consturct the actual P1 function here
  %  nodeValuesJ1 should be the same as vJ1(n4e), so reverse that
end
