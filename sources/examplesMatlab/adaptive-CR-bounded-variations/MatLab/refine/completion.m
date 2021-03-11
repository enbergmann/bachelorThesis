function [c4n,n4e,n4sDb,n4sNb,err4e]=completion(c4n_nc,...
    c4n, n4e, n4sDb,n4sNb,err4e,rhsf)
%% completion - generate a regular triangulation
%               by means of completion algorithm
%               needed for the real tsa + completion algorithm
%
%% Input
% c4n_nc,n4e_nc,n4sDb_nc,n4sNb_nc 
%               - triangulation with possible hanging nodes
% c4n, n4e, n4sDb, n4sNB          
%               - last regular triangulation
% err4e,rhsf    - input to compute e and etilde of tsa
%
%% Output
% c4n,n4e,n4sDb,n4sNb - regular triangulation after completion
% err4e               - refinement indicator for tsa (e, etilde)

%    Copyright (C) 2013  Hella Rabus
%                   
%    You should have received LICENSE.txt along with this file
%    that gives further information an the license.
%    See the GNU General Public License for more details.


% compute all midpoints the last regular triangulation triangulation
n4s = computeN4s(n4e);
mid4s = computeMid4s(c4n, n4s);
% index = numbers the midpoints, that belong to the _nc triangulation 
% (i.e. these edges will be refined)

% Mark
c4n_new = unique(c4n_nc,'rows');
c4n_new = setdiff(c4n_new, c4n, 'rows');
[~, index, indexc4n] = intersect(mid4s, c4n_new, 'rows');
fprintf(1,'     completion: nr Nodes in irregular mesh: %6.0f \n',...
    size(unique(c4n_nc, 'rows'), 1));
fprintf(1,'     nr Nodes: %9.0f; nr of marked edges: %7.0f',...
    size(c4n, 1),size(index, 1));
while ~isempty(index)&&~isempty(c4n_new)

    % REFINE
    [c4n, n4e, n4sDb, n4sNb, err4e] = refineBi3GB_irregular(c4n, n4e,...
         n4sDb, n4sNb, n4s(index, :), err4e, rhsf, true);
    c4n_new = setdiff(c4n_new, c4n_new(indexc4n, :), 'rows');    
    %MARK
    n4s = computeN4s(n4e);
    mid4s = computeMid4s(c4n, n4s);
    % Update the Output on the screen
    [~,index,indexc4n] = intersect(mid4s,c4n_new,'rows');
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf(1,'\b');
    fprintf(1,'nr Nodes: %9.0f; nr of marked edges: %7.0f\n',...
        size(c4n, 1), size(index, 1));
end
