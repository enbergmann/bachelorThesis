function absNodeJumps4s = computeAbsNodeJumps4s(n4e, e4s, nodeValues4e)
% Computes the jumps of a elementwise defined function with local nodal values
% given by nodeValues4e in the two nodes of each side of the triangulation.
%
% computeAbsNodeJumps4s.m
% input:  n4e            - nodes for elements
%         e4s            - elements for sides 
%         nodeValues4e   - '(size(n4e))-dimensional double array' where the
%                          k-th entry of the j-th row contains the value in the
%                          k-th local node of the j-th triangle of the mesh
%
% output: absNodeJumps4s - '(size(e4s))-dimensional double array' where the
%                          j-th row contains the two jumps in the two nodes of
%                          the j-th side
  absNodeJumps4s = zeros(size(e4s));
  for side = 1:size(e4s, 1)
    tPlus = e4s(side, 1);
    tMinus = e4s(side, 2);
    if tMinus ~= 0 
      % if there are two triangles with side as an edge, i.e.  side is an inner
      % edge (else there are no jumps)
      pos = 1;
      for indPlus = 1:3
        indMinus = find(n4e(tMinus, :) == n4e(tPlus, indPlus));
          % is the indPlus-th node of tPlus a node of tMinus and which one    
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
