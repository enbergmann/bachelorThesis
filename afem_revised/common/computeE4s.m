function e4s = computeE4s(n4e)
%% computeE4s - Elements for sides.
%   computeE4s(n4e) returns a matrix each row of which corresponds to one side
%               of the decomposition. The side numbering is the same as in
%               n4s. Each row contains the numbers of the two elements that
%               the corresponding side is a part of. If it is a boundary
%               side the second entry is 0.
%               n4e is as specified in the documentation.
%
%   See also: computeN4s
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
% Modified (P. Bringmann): local numbering of edges s.t. P_k is opposite
% node of E_k

    if isempty(n4e)
        e4s = zeros(0,2);
        return;
    end

    %% Compute e4s.
    allSides = [n4e(:,[2 3]); n4e(:,[3 1]); n4e(:,[1 2])];
    [b,ind,back] = unique(sort(allSides,2),'rows','first');
    n4sInd = sort(ind); % by the way: n4s = allSides(n4sInd,:)

    nrElems = size(n4e,1);
    elemNumbers = [1:nrElems 1:nrElems 1:nrElems];
    e4s(:,1) = elemNumbers(n4sInd);
    allElem4s(ind) = accumarray(back,elemNumbers);
    e4s(:,2) = allElem4s(n4sInd)'-e4s(:,1);
end
