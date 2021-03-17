function [c4n, n4e, n4sDb, n4sNb] = overlay(c4n, n4e, n4sDb, ...
            n4sNb,c4n2,c4n3)

%% overlay - Computes the overlay (coarsest common refinement) 
%            of two triangulations T2 and T3 refined by NVB
%            from a coarse triangulation T0
% INPUT 
%     c4n, n4e, n4sDb, n4sNb 
%           - initial (regular) triangulation T0
%       c4n2     - (regular) triangulation T2 refined from T0
%       c4n3     - (regular)triangulation T3 refined from T0
%
% OUTPUT
%     c4n, n4e, n4sDb, n4sNb
%                - overlay triangulation of T2 and T3

%    Copyright (C) 2013  Hella Rabus
%                   
%    You should have received LICENSE.txt along with this file
%    that gives further information an the license.
%    See the GNU General Public License for more details.


% matrix of coordinates of all nodes that will belong 
% to the overlay triangulation
c = unique([c4n2;c4n3],'rows');

%% MARK
% compute n4s and all midpoints in T0; 
n4s = computeN4s(n4e);
mid4s = computeMid4s(c4n,n4s);

% index = numbers the midpoints of edges in T0,
%    that belong to the overlay 
%    (i.e. these edges will be refined)
[dummy{2},index] = intersect(mid4s,c,'rows');

% output
fprintf(1,'\n compute overlay:');

% Loop: 
% as long as there are coordinates in c that match a midpoint
while ~isempty(index)
    
    % output points, that indicate further computing
    fprintf(1,'.');
    
    % REFINE the edges, whose midpoints belong to the overlay
    [c4n, n4e, n4sDb, n4sNb] = refineBi3GB(c4n, n4e, n4sDb, ...
            n4sNb,n4s(index,:));
    %compute n4s and all midpoints
    n4s = computeN4s(n4e);
    mid4s = computeMid4s(c4n,n4s);
    % MARK : the intersection of c and the midpoints defines 
    %    the set of to be refined edges
    [dummy{2},index] = intersect(mid4s,c,'rows');
end
fprintf(1,'\n');
