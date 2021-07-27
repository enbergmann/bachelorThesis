function n4sMarked = closure(n4e, n4sMarked)
%% closure - Mark Reference Sides.
%   closure(n4e, n4sMarked) markes the reference side of each element with
%       a marked side. The reference side of an element is its first one.
%   	n4e is as specified in the documentation, n4sMarked has the same
%       structure as n4s.
%       The output is a new list of marked sides having the same structure
%       as n4s.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

  %% INITIALIZATION
  n4s = computeN4s(n4e);
  s4n = computeS4n(n4e, n4s);
  s4e = computeS4e(n4e);
  nSides = size(n4s, 1);

  MarkedSides =...
    full(s4n(sub2ind(size(s4n), n4sMarked(:,1), n4sMarked(:,2))));

  isMarked = zeros(nSides, 1);
  isMarked(MarkedSides) = 1;
  isMarked4e = isMarked(s4e);

  %% CLOSURE
  while true
    % Check for consistency
    % isConsistent(k) == 0 if element k has a marked side, but its
    %                      reference side is not marked
    %                 == 1 otherwise
    isConsistent = ~and(any(isMarked4e(:,[1 2]),2), ~isMarked4e(:,3));
    if all(isConsistent)  % if all elements are consistent => done!
        break;
    end
    % Mark sides that need to be refined
    n4sMarked_new = n4e(~isConsistent,[1 2]);
    % Remove duplettes
    n4sMarked_new = unique(sort(n4sMarked_new, 2), 'rows', 'stable');
    n4sMarked = [n4sMarked; n4sMarked_new];

    NewMarkedSides =...
      full(s4n(sub2ind(size(s4n), n4sMarked_new(:,1), n4sMarked_new(:,2))));

    isMarked(NewMarkedSides) = 1;
    isMarked4e = isMarked(s4e);

  end
end

