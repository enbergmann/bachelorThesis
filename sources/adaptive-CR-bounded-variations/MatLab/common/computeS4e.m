function s4e = computeS4e(n4e)
%% computeS4e - Sides for elements.
%   computeS4e(n4e) returns a matrix each row of which corresponds to one
%               element of the decomposition. Each row contains the numbers
%               of the three sides belonging to an element. The side
%               numbering is the same as in n4s.
%
%   See also: computeN4s
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    if isempty(n4e)
        s4e = [];
        return;
    end

    %% Compute s4e.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    [b,ind,back] = unique(sort(allSides,2),'rows','first');
    [n4sInd, sortInd] = sort(ind); % by the way: n4s = allSides(n4sInd,:)
    sideNr(sortInd) = 1:length(ind); % sideNr(back): numbers for allSides
    s4e = reshape(sideNr(back),size(n4e));
end