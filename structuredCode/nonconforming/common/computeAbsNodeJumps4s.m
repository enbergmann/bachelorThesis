function absNodeJumps4s = computeAbsNodeJumps4s(n4e, e4s, nodeValues4e)
%% DOC
% Computes the jumps of a elementwise defined function u with local nodal
% values given by nodeValues4e in the two nodes of each side of the
% triangulation.
%
% computeAbsNodeJumps4s.m
% input: n4e          - nodes for elements
%        e4s          - elements for sides 
%        nodeValues4e - '(nrElems x 3)-dimensional double array' where the k-th
%                       entry of the j-th row contains the value of u in the
%                       k-th local node of the j-th triangle of the mesh
%
% output: absNodeJumps4s - '(nrElems x 2)-dimensional double array' where the
%                          j-th row contains the two jumps in the two nodes of
%                          the j-th side of the triangulation

%% INIT
  absNodeJumps4s = zeros(size(e4s));

%% MAIN
  for side = 1:size(e4s, 1)
    tPlus = e4s(side, 1);
    tMinus = e4s(side, 2);
      % the two neighbouring triangles of side (tMinus = 0 if side is an outer
      % edge)
    if tMinus ~= 0 
      pos = 1;
      for indPlus = 1:3
        indMinus = find(n4e(tMinus, :) == n4e(tPlus, indPlus));
          % find node of tMinus that is the indPlus-th node of tPlus (if it is
          % one at all)
        if not(isempty(indMinus))
          % if the indPlus-th node of tPlus is the indMinus-th node of
          % tMinus...
          absNodeJumps4s(side, pos) = abs(nodeValues4e(tPlus, indPlus) - ...
                                         nodeValues4e(tMinus, indMinus));
            % ... compute the jump in the pos-th node of side
          pos = pos + 1;
        end
      end
    end
  end
end
