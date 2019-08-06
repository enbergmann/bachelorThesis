function val = rowaddr(A,I,J)
%% rowaddr
% rowaddr(A,I,J) returns indices val for Matrix A and subsets I of
% A(:,1) and J of (A:,2) such that A(val(k),:) = [I(k),J(k)]
%
% (C) 2009--2013 C. Merdon, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

val = full(A(sub2ind(size(A),I(:),J(:))));
end
