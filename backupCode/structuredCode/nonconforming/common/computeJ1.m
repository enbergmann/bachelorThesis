function vJ1 = computeJ1(n4e, n4sB, v)
%% DOC
% Computes the enriching operator of a CR function v defined by 
%   J_1 v (z) := |\Tcal(z)|^{-1} \sum_{T\in\Tcal(z)}v|_T(z) 
% for all z\in\Ncal(\Omega) as a P_1(\Tcal)\cap C_0(\overline\Omega) (Courant)
% function with respect to the mesh defined by [c4n, n4e].
%
% computeJ1.m
% input: c4n  - coordinates for nodes
%        n4e  - nodes for elements
%        n4sB - nodes for boundary sides
%        v    - '(nrSides x 1)-dimensional double array' where the j-th row
%               contains the coefficient of the CR function w.r.t. the j-th
%               side of the triangulation
%
% output: vJ1 - '(nrNodes x 1)-dimensional double array' where the j-th row
%               contains the coefficient of J1 v w.r.t. the j-th node of the
%               triangulation

%% INIT
  s4e = computeS4e(n4e); 
  nodeValuesCR4e = computeNodeValuesCR4e(s4e, v);

  nrNodes = max(max(n4e));
  innerNodes = setdiff(1:nrNodes, unique(n4sB(:)));

%% MAIN
  % compute enriching operator
  nodeValuesJ1 = zeros(size(nodeValuesCR4e));
  for j = innerNodes
    indNode = find(n4e==j);
    localVals4node = nodeValuesCR4e(indNode);
    nodeValuesJ1(indNode) = sum(localVals4node)/length(localVals4node);
  end

  vJ1 = zeros(nrNodes, 1);
  vJ1(n4e) = nodeValuesJ1;
end
