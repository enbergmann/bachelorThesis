function absNodeJumps4s = ...
    computeAbsNodeJumps4s(n4e, n4s, e4s, nrSides, nodeValues4e)
%% DOC
% Computes the absolute jumps of a element-wise defined function v with local
% nodal values given by nodeValues4e in the two nodes of each side of the
% triangulation.
%
% computeAbsNodeJumps4s.m
% input: n4e          - nodes for elements
%        n4s          - nodes for sides
%        e4s          - elements for sides 
%        nrSides      - number of sides
%        nodeValues4e - '(nrElems x 3)-dimensional double array' where the k-th
%                       entry of the j-th row contains the value of v in the
%                       k-th local node of the j-th triangle of the mesh
%
% output: absNodeJumps4s - '(nrElems x 2)-dimensional double array' where the
%                          j-th row contains the two absolute jumps in the two
%                          nodes of the j-th side of the triangulation (if a
%                          side is an outer edge the jump in a node is the 
%                          value of v in that node)

%% INIT
  absNodeJumps4s = zeros(size(e4s));

%% MAIN
  for side = 1:nrSides
    tPlus = e4s(side, 1);
    tMinus = e4s(side, 2);
      % the two neighbouring triangles of side (tMinus = 0 if side is an outer
      % edge)
    node1 = n4s(side, 1);
    node2 = n4s(side, 2);
      % the two nodes of side = conv({node1, node2})
    indPlus1 = n4e(tPlus, :) == node1;
    indPlus2 = n4e(tPlus, :) == node2;
      % find the local node numbers of node1 and node2 on tPlus
    valPlus1 = nodeValues4e(tPlus, indPlus1);
    valPlus2 = nodeValues4e(tPlus, indPlus2);
      % the values of v in node1 and node2 on tPlus
    if tMinus == 0
      valMinus1 = 0;
      valMinus2 = 0;
        % the values of v in node1 and node2 on tMinus if side is an outer edge
    else
      indMinus1 = n4e(tMinus, :) == node1;
      indMinus2 = n4e(tMinus, :) == node2;
        % find the local node numbers of node1 and node2 on tMinus if side
        % is an inner edge
      valMinus1 = nodeValues4e(tMinus, indMinus1);
      valMinus2 = nodeValues4e(tMinus, indMinus2);
        % the values of v in node1 and node2 on tMinus if side is an inner edge
    end
    absNodeJumps4s(side, :) = abs([valPlus1-valMinus1, valPlus2 - valMinus2]);
  end
end
