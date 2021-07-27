function n4s = refineEdges(n4s, n4sMarked, s4n)
% Refine the edges by splitting each marked side into two.

% (C) 2016 Bringmann

  if isempty(n4s)
    return
  end

  %% INITIALIZATION
  nNodes = size(s4n, 1);
  nSides = max(nonzeros(s4n));

  Sides = full(s4n(sub2ind(size(s4n), n4s(:,1), n4s(:,2))));
  MarkedSides =...
    full(s4n(sub2ind(size(s4n), n4sMarked(:,1), n4sMarked(:,2))));
  nMarkedSides = length(MarkedSides);

  % newNodes4s(k) == j > 0 if side k is marked with new node j
  %               == 0     otherwise
  newNodes4s = zeros(nSides, 1);
  newNodes4s(MarkedSides) = nNodes + (1:nMarkedSides);
  newNodes4s = newNodes4s(Sides);

  %% MERGE
  n4s = [n4s(:,1),   newNodes4s, newNodes4s, n4s(:,2)]';
  n4s = reshape(n4s(n4s ~= 0), 2, [])';

end

