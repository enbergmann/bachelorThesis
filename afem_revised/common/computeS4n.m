function s4n = computeS4n(n4e,n4s)
%% computeS4n
% computeS4n(n4e[,n4s]) returns a symmetric sparse matrix in
% which the entry (j,k) contains the number of the side with
% the end nodes j and k or zero if no such side exists.
% The side enumeration is the same as in n4s.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

  if nargin<2
      n4s = computeN4s(n4e);
  end
  S = size(n4s,1);
  N = max(n4e(:));
  s4n = sparse(n4s(:,1),n4s(:,2),1:S,N,N);
  % Up to here, s4n is not yet symmetric as each side has only
  % been considered once. The following makes sure that s4n is
  % symmetric.
  s4n = s4n + s4n';
end