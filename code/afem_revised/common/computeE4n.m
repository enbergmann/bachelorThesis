function e4n = computeE4n(n4e)
%% computeE4n - Elements for nodes.
%   computeE4n(n4e) returns a sparse matrix in which the entry (j,k) contains
%               the number of the element that has the nodes j and k in
%               counterclockwise order (or 0 if no such element exists).
%               n4e is as specified in the documentation.
%
%   See also: computeE4s
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    if isempty(n4e)
        e4n = [];
        return;
    end

    %% Compute e4n.
    % Create a list of all sides in the decomposition and build a sparse
    % matrix such that each side computes its proper element number.
    allSides = [n4e(:,[2 3]); n4e(:,[3 1]); n4e(:,[1 2])];
    nrElems = size(n4e,1);
    N = max(max(n4e));
    elemNumbers = [1:nrElems 1:nrElems 1:nrElems];
    e4n = sparse(allSides(:,1),allSides(:,2),elemNumbers,N,N);
end
