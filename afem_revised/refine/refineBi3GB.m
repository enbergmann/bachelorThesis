function [c4n, n4e, n4sDb, n4sNb, parents4e, child4e,...
          nChildren4e, childNo4e] =...
  refineBi3GB(c4n, n4e, n4sDb, n4sNb, n4sMarked, allowHangingNodes)
%% refineBi3GB - refine using the Bisec3-Green-Blue-strategy
%   refineBi3GB(c4n, n4e, n4sDb, n4sNb, n4sMarked) Refines a given mesh using
%       the Bisec3-Green-Blue refinement. For details on data structures
%       and refinement strategies see the documentation. Input is a mesh
%       defined by c4n, n4e, n4sDb, n4sNb and marked sides given by
%       n4sMarked. Output is a refined mesh defined by c4nNew, n4eNew,
%       n4sDbNew and n4sNbNew.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

  %% PROCEED INPUT
  if nargin < 6
    allowHangingNodes = false;
  end

  %% INITIALIZATION
  n4s = computeN4s(n4e);
  s4n = computeS4n(n4e, n4s);
  s4e = computeS4e(n4e);
  nElem = size(n4e, 1);
  nSides = size(n4s, 1);
  nNodes = size(c4n, 1);

  %% CLOSURE
  % n4sMarked = closure(n4e, n4sMarked);

  %% UPDATE C4N
  c4n = [c4n; computeMid4s(c4n, n4sMarked)];

  %% COMPUTE EDGES TO REFINE
  if allowHangingNodes
      % new nodes numbered according to n4sMarked
    [side1Marked4e, side1NewNodeNo4e] = ...
      ismember(n4e(:,[2 3]), n4sMarked, 'rows');
    [side2Marked4e, side2NewNodeNo4e] = ...
      ismember(n4e(:,[3 1]), n4sMarked, 'rows');
    [side3Marked4e, side3NewNodeNo4e] = ...
      ismember(n4e(:,[1 2]), n4sMarked, 'rows');
    isMarked4e = [side1Marked4e, side2Marked4e, side3Marked4e];
    newNodes4e = zeros(nElem, 3);
    newNodes4e(isMarked4e) = nNodes + ...
                             nonzeros([side1NewNodeNo4e;
                                       side2NewNodeNo4e;
                                       side3NewNodeNo4e]);
  else
    MarkedSides =...
      full(s4n(sub2ind(size(s4n), n4sMarked(:,1), n4sMarked(:,2))));

    % newNodes4s(k) == j > 0 if side k is marked with new node j
    %               == 0     otherwise
    newNodes4s = zeros(nSides, 1);
    newNodes4s(MarkedSides) = nNodes + (1:length(MarkedSides));

    % newNodes4e(k,m) == j > 0 if side m (1,2 or 3) of (old) element k is
    %                          marked with new node j
    %                 == 0     if side m of (old) element k is not marked
    newNodes4e = newNodes4s(s4e);
  end

  %% REFINEMENT ALGORITHM
  %        n3        Element [n1 n2 n3] (as of a row in n4e)
  %       /  \       has the sides/new nodes [s1 s2 s3] (as
  %     s2    s1     given by a row in newNodes4e), in the order
  %    /        \    depicted here.  Each sj can be a new
  %   n1 --s3-- n2   node, depending of the refinement type.

  % Get lists of elements for corresponding refinement rules
  % Not to be refined
  Ind0 = ~any(newNodes4e,2);
  % To be bisec3 refined (all sides are marked)
  IndBi3 = all(newNodes4e,2);
  % To be blue_right refined (reference and 1st edge is marked)
  IndB1 = and(all(newNodes4e(:,[1 3]),2),~newNodes4e(:,2));
  % To be blue_left refined (reference and 2nd edge is marked)
  IndB2 = and(all(newNodes4e(:,[2 3]),2),~newNodes4e(:,1));
  % To be green refined (only reference edge is marked)
  IndG = and(newNodes4e(:,3),~any(newNodes4e(:,[1 2]),2));

  % Remaining cases for irregular refinement
  % Irregular refinement of edge 1 (green)
  Irr1 = and(newNodes4e(:,1),~any(newNodes4e(:,[2 3]),2));
  % Irregular refinement of edge 2 (green)
  Irr2 = and(newNodes4e(:,2),~any(newNodes4e(:,[1 3]),2));
  % Irregular refinement of edges 1 and 2 (blue right)
  Irr12 = and(all(newNodes4e(:,[1 2]),2),~newNodes4e(:,3));

  if any([Irr1; Irr2; Irr12])
    warning('Irregular refinement in refineBi3GB')
  end

  n4e = [% Untouched
         n4e(Ind0,:)
         % Bisec3
         %       o
         %      /:\
         %     / : \
         %    o 2:3 o
         %   / \ : / \
         %  / 1 \:/ 4 \
         % o-----o-----o
         n4e(IndBi3,1)         newNodes4e(IndBi3,3)  newNodes4e(IndBi3,2)
         newNodes4e(IndBi3,3)  n4e(IndBi3,3)         newNodes4e(IndBi3,2)
         n4e(IndBi3,3)         newNodes4e(IndBi3,3)  newNodes4e(IndBi3,1)
         newNodes4e(IndBi3,3)  n4e(IndBi3,2)         newNodes4e(IndBi3,1)
         % Blue Right
         %       o
         %      /:\
         %     / : \
         %    /  :3 o
         %   / 1 : / \
         %  /    :/ 2 \
         % o-----o-----o
         n4e(IndB1,3)        n4e(IndB1,1)        newNodes4e(IndB1,3)
         newNodes4e(IndB1,3) n4e(IndB1,2)        newNodes4e(IndB1,1)
         n4e(IndB1,3)        newNodes4e(IndB1,3) newNodes4e(IndB1,1)
         % Blue Left
         %       o
         %      /:\
         %     / : \
         %    o 2:  \
         %   / \ : 3 \
         %  / 1 \:    \
         % o-----o-----o
         n4e(IndB2,1)        newNodes4e(IndB2,3) newNodes4e(IndB2,2)
         newNodes4e(IndB2,3) n4e(IndB2,3)        newNodes4e(IndB2,2)
         n4e(IndB2,2)        n4e(IndB2,3)        newNodes4e(IndB2,3)
         % Green
         %       o
         %      /:\
         %     / : \
         %    /  :  \
         %   / 1 : 2 \
         %  /    :    \
         % o-----o-----o
         n4e(IndG,3)         n4e(IndG,1)         newNodes4e(IndG,3)
         n4e(IndG,2)         n4e(IndG,3)         newNodes4e(IndG,3)
         % Irregular refinement of side 1
         %       o
         %      / \
         %     /   \
         %    / 2  .o
         %   /  .'   \
         %  /.'   1   \
         % o-----------o
         n4e(Irr1,1)         n4e(Irr1,2)         newNodes4e(Irr1,1)
         n4e(Irr1,3)         n4e(Irr1,1)         newNodes4e(Irr1,1)
         % Irregular refinement of side 2
         %       o
         %      / \
         %     /   \
         %    o.  2 \
         %   /  ' .  \
         %  /  1   ' .\
         % o-----------o
         n4e(Irr2,1)         n4e(Irr2,2)         newNodes4e(Irr2,2)
         n4e(Irr2,3)         n4e(Irr2,1)         newNodes4e(Irr2,2)
         % Irregular refinement of side 1 and 2
         %       o
         %      / \
         %     / 2 \
         %    o----.o
         %   / 3.'   \
         %  /.'   1   \
         % o-----------o
         n4e(Irr12,1)        n4e(Irr12,2)        newNodes4e(Irr12,1)
         newNodes4e(Irr12,1) n4e(Irr12,3)        newNodes4e(Irr12,2)
         n4e(Irr12,1)        newNodes4e(Irr12,1) newNodes4e(Irr12,2)
        ];

  %% REFINE BOUNDARY
  n4sDb = refineEdges(n4sDb, n4sMarked, s4n);
  n4sNb = refineEdges(n4sNb, n4sMarked, s4n);


  %% RECOVERING STRUCTURE FOR OLD DATA
  if nargout > 4
    Elements = 1:nElem;
    parents4e = [repmat(Elements(Ind0),     1, 1)
                 repmat(Elements(IndBi3),   4, 1)
                 repmat(Elements(IndB1),    3, 1)
                 repmat(Elements(IndB2),    3, 1)
                 repmat(Elements(IndG),     2, 1)
                 repmat(Elements(IndIrr1),  2, 1)
                 repmat(Elements(IndIrr2),  2, 1)
                 repmat(Elements(IndIrr12), 3, 1)
                 ];
  end

  if nargout > 5
    nElem0 = nnz(Ind0);
    nElemBi3 = nnz(IndBi3);
    nElemB1 = nnz(IndB1);
    nElemB2 = nnz(IndB2);
    nElemG = nnz(IndG);
    nElemIrr1 = nnz(IndIrr1);
    nElemIrr2 = nnz(IndIrr2);
    nElemIrr12 = nnz(IndIrr12);

    % TODO explain duplettes
    child4e = [ ...
       % unrefined
       repmat((1:nElem0)', 1, 4);...
       % bisec3
       (1:nElemBi3)'   + nElem0,...
       (1:nElemBi3)'   + nElem0 +   nElemBi3,...
       (1:nElemBi3)'   + nElem0 + 2*nElemBi3,...
       (1:nElemBi3)'   + nElem0 + 3*nElemBi3;...
       % blue_right
       (1:nElemB1)'    + nElem0 + 4*nElemBi3,...
       (1:nElemB1)'    + nElem0 + 4*nElemBi3,...
       (1:nElemB1)'    + nElem0 + 4*nElemBi3 +   nElemB1,...
       (1:nElemB1)'    + nElem0 + 4*nElemBi3 + 2*nElemB1;...
       % blue_left
       (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1,...
       (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1 +   nElemB2,...
       (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1 + 2*nElemB2,...
       (1:nElemB2)'    + nElem0 + 4*nElemBi3 + 3*nElemB1 + 2*nElemB2;...
       % green
       (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2,...
       (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2,...
       (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       +   nElemG,...
       (1:nElemG)'     + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       +   nElemG;...
       % irregular refinement of edge 1
       (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG,...
       (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG,...
       (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG +   nElemIrr1,...
       (1:nElemIrr1)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG +   nElemIrr1;...
       % irregular refinement of edge 2
       (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1,...
       (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1,...
       (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1 +   nElemIrr2,...
       (1:nElemIrr2)'  + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1 +   nElemIrr2;...
       % irregular refinement of edge 1 and 2
       (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2,...
       (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2,...
       (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2 + nElemIrr12,...
       (1:nElemIrr12)' + nElem0 + 4*nElemBi3 + 3*nElemB1 + 3*nElemB2...
                       + 2*nElemG + 2*nElemIrr1 + 2*nElemIrr2 + nElemIrr12
      ];
  end

  if nargout > 6
    nChildren4e = [  ones(nElem0,     1)
                   4*ones(nElemBi3,   1)
                   3*ones(nElemB1,    1)
                   3*ones(nElemB2,    1)
                   2*ones(nElemG,     1)
                   2*ones(nElemIrr1,  1)
                   2*ones(nElemIrr2,  1)
                   3*ones(nElemIrr12, 1)
                   ];
  end

  if nargout > 7
    childNo4e = [  ones(nElem0,     1)
                   ones(nElemBi3,   1)
                 2*ones(nElemBi3,   1)
                 3*ones(nElemBi3,   1)
                 4*ones(nElemBi3,   1)
                   ones(nElemB1,    1)
                 2*ones(nElemB1,    1)
                 3*ones(nElemB1,    1)
                   ones(nElemB2,    1)
                 2*ones(nElemB2,    1)
                 3*ones(nElemB2,    1)
                   ones(nElemG,     1)
                 2*ones(nElemG,     1)
                   ones(nElemIrr1,  1)
                 2*ones(nElemIrr1,  1)
                   ones(nElemIrr2,  1)
                 2*ones(nElemIrr2,  1)
                   ones(nElemIrr12, 1)
                 2*ones(nElemIrr12, 1)
                 3*ones(nElemIrr12, 1)
                 ];
  end

end

