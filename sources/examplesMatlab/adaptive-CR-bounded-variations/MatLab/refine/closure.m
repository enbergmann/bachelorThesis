function n4sRefine = closure(n4e, n4sMarked)
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

    %% Preliminary work.
    nrNodes = max(max(n4e));
    nrElems = size(n4e,1);
    e4n = computeE4n(n4e);

    %% Initialize the matrix of marked sides.
    %  This is a symmetric matrix where entry (j,k) is 1 iff there is a
    %  marked side between nodes j and k.
    marked4n = sparse(n4sMarked(:,[1 2]), n4sMarked(:,[2 1]),...
                      ones(size(n4sMarked,1),2),nrNodes,nrNodes);

    %% Execute the closure algorithm.
    for curElem = 1 : nrElems
        curNodes = n4e(curElem,:);
        % While the current element's reference side is not marked and
        % another side of the current element is marked....
        while marked4n(curNodes(1),curNodes(2)) == 0 ...
              && ( marked4n(curNodes(2),curNodes(3)) == 1 ...
                   || marked4n(curNodes(3),curNodes(1)) == 1 )
            % ... mark the reference side.
            marked4n(curNodes(1),curNodes(2)) = 1;
            marked4n(curNodes(2),curNodes(1)) = 1;
            % Marking the reference side of the current element may have
            % created a similar situation on the neighbouring element. We
            % need to handle this.
            % Tricky but true: e4n(curNodes(1),curNodes(2)) will return
            % curElem while this returns the neighbouring element number:
            neighbourElem = e4n(curNodes(2),curNodes(1));
            if neighbourElem > 0
                curNodes = n4e(neighbourElem,:);
            end
        end
    end

    %% Assemble the side list for refinement.
    marked4n = triu(marked4n);
    [row,col] = find(marked4n);
    n4sRefine = [row,col];
end

