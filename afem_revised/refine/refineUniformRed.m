function [c4n,n4e,n4sDb,n4sNb,parents4e,pos4pelems,parents4Dbs,parents4Nbs,childnr]...
    = refineUniformRed(c4n,n4e,n4sDb,n4sNb)
%% refineUniformRed
% red mesh refinement by connection of side midpoints
% (--> 4 subtriangles with same area)
%
% Input:     c4n          coordinates for the nodes of the mesh
%            n4e          nodes for the elements of the mesh
%            n4sDb        the nodes of the sides in the Dirichlet boundary
%            n4sNb        the nodes of the sides in the Neumann boundary
%
% Output:    c4n          ---\
%            n4e              data of red-refined triangulation
%            n4sDb        ---/
%            n4sNb        --/
%            parents4e    parent triangle for fine triangle
%            parents4sDbs parent edge for fine Dirichlet edges
%            parents4sNbs parent edge for fine Neumann edges
%            childnr      child number type (1,2,3 or 4) for fine elements
%
% mostly (C) 2009--2012 W. Boiger, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
% Modified (C. Merdon): treatment of boundary data, some fixes and comments
% Modified (P. Bringmann): local numbering of edges s.t. P_k is opposite
% node of E_k

n4s = computeN4s(n4e);
s4n = computeS4n(n4e);
s4e = computeS4e(n4e);
mid4s = computeMid4s(c4n, n4s);
nNodes = size(c4n, 1);
nElem = size(n4e, 1);
% Add new nodes to c4n
c4n = [c4n; mid4s];
clear('mid4s');
% Store relationship between old and new nodes
parents4e = [1:nElem 1:nElem 1:nElem 1:nElem];
clear('n4s');
% Core routine: Create new elements
%      n3       Element [n1 n2 n3] (as of a row in n4e)
%     /  \      has the sides [s1 s2 s3] (as given by a
%   s2    s1    row in s4e), in the order depicted here.
%  /       \    Each side corresponds to a new node ---
% n1 --s3-- n2  this node is the midpoint of the side.
s4e = s4e + nNodes;
n4e = [n4e(:,1) s4e(:,3) s4e(:,2);  % [n1 s3 s2]
       s4e(:,3) n4e(:,2) s4e(:,1);  % [s3 n2 s1]
       s4e(:,2) s4e(:,1) n4e(:,3);  % [s2 s1 n3]
       s4e];                        % [s1 s2 s3]
childnr = [  ones(nElem, 1);
           2*ones(nElem, 1);
           3*ones(nElem, 1);
           4*ones(nElem, 1)];
% child numbers
%       o
%      / \
%     / 3 \
%    o-----o
%   / \ 4 / \
%  / 1 \ / 2 \
% o-----o-----o
pos4pelems = reshape([1 2 3 0]'*ones(1, nElem), [1 4*nElem]);
clear('s4e');
% Refine boundary
parents4Dbs = reshape([1:size(n4sDb,1); 1:size(n4sDb)], [1 2*size(n4sDb,1)]);
parents4Nbs = reshape([1:size(n4sNb,1); 1:size(n4sNb)], [1 2*size(n4sNb,1)]);
n4sDb = refineBoundary(n4sDb, s4n, nNodes);
n4sNb = refineBoundary(n4sNb, s4n, nNodes);
end

function n4sB = refineBoundary(n4sB, s4n, offset)
% Refine the boundary by splitting each side into two.
if isempty(n4sB)
    return
end
% Get new (middle) node for each boundary side
midB = full(s4n(sub2ind(size(s4n), n4sB(:,1), n4sB(:,2))));
midB = midB + offset;
% Split every side by inserting the new nodes (twice)
midB = reshape(midB, 1, []);
n4sB = n4sB';
n4sB = [n4sB(1,:); midB; midB; n4sB(2,:)];
n4sB = reshape(n4sB, 2, []);
n4sB = n4sB';
end
