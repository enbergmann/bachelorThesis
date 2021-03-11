function [c4nNew,n4eNew,n4sDbNew,n4sNbNew] = refineUniformRed(c4n,n4e,n4sDb,n4sNb)
%% refineUniformRed - Refine every element the "red" way.
%   refineUniformRed(c4n, n4e, n4sDb, n4sNb) Refines a given mesh uniformly
%       using the red refinement. For details on data structures and
%       refinement strategies see the documentation. Input is a mesh
%       defined by c4n, n4e, n4sDb and n4sNb. Output is a refined mesh
%       defined by c4nNew, n4eNew, n4sDbNew and n4sNbNew.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    %% Preliminary work.
    nrNodes = size(c4n,1);
    nrElems = size(n4e,1);
    n4s = computeN4s(n4e);
    nrSides = size(n4s,1);
    newNodes4s = sparse(n4s(:,1),n4s(:,2),(1:nrSides)'+ nrNodes, ...
                        nrNodes,nrNodes);
    newNodes4s = newNodes4s + newNodes4s';

    %% Compute coordinates for new nodes.
    % As every element is refined "red", there is a new node on each side.
    mid4s = computeMid4s(c4n,n4s);
    c4nNew = [c4n;mid4s];

    %% red refinement
    n4eNew = zeros(4*nrElems,3);
    for curElem = 1 : nrElems
        curNodes = n4e(curElem,:);
        curNewNodes = [newNodes4s(curNodes(1),curNodes(2));
                       newNodes4s(curNodes(2),curNodes(3));
                       newNodes4s(curNodes(3),curNodes(1));
                      ];
        % Generate new elements.
        n4eNew(4*(curElem-1)+1:4*curElem,:) = ...
                  [ curNodes(1)    curNewNodes(1) curNewNodes(3);
                    curNewNodes(1) curNodes(2)    curNewNodes(2);
                    curNewNodes(2) curNewNodes(3) curNewNodes(1);
                    curNewNodes(3) curNewNodes(2) curNodes(3);
                  ];
    end

    %% refinement of  Dirichlet boundary
    n4sDbNew = zeros(2*size(n4sDb,1),2);
    for curSide = 1 : size(n4sDb,1)
        curNodes = n4sDb(curSide,:);
        curNewNodes = newNodes4s(curNodes(1),curNodes(2));
        % Generate new Dirichlet boundary sides.
        n4sDbNew(2*(curSide-1)+1:2*curSide,:) = ...
                  [ curNodes(1)    curNewNodes;
                    curNewNodes    curNodes(2);
                  ];
    end

    %% refinement of  Neumann boundary
    n4sNbNew = zeros(2*size(n4sNb,1),2);
    for curSide = 1 : size(n4sNb,1)
        curNodes = n4sNb(curSide,:);
        curNewNodes = newNodes4s(curNodes(1),curNodes(2));
        % Generate new Neumann boundary sides.
        n4sNbNew(2*(curSide-1)+1:2*curSide,:) = ...
                  [ curNodes(1)    curNewNodes;
                    curNewNodes    curNodes(2);
                  ];
    end
end
