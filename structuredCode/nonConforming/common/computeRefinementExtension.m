% TODO comment and structure and write the interface documentation
% function computeRefinementExtension
%
% Compute the prolongation of a solution of the CRFEM to a
% one-level refinement of the mesh.
%
% INPUT         c4n          Coordinates of the coarse mesh
%               n4e          nodes for elements of the coarse mesh
%               c4nNew          Coordinates of the fine mesh
%               n4eNew          nodes for elements of the fine mesh
%               u          solution on the coarse mesh
%
% OUTPUT        uNew          prolongation of the old solution to the new
%                               mesh

function uNew = computeRefinementExtension(c4n, n4e, c4nNew, n4eNew, u)
  %% initialization
  n4sNew = computeN4s(n4eNew);
  s4n = computeS4n(n4e);
  s4e = computeS4e(n4e); 
  s4eNew = computeS4e(n4eNew);

  nrNodes = size(c4n, 1); 
  nrNodesNew = size(c4nNew, 1);
  nrElemsNew = size(n4eNew, 1);
  nrSidesNew = size(n4sNew, 1);

  n4parentSides4n = getParentSide(c4n, c4nNew, s4n, n4sNew, ...
    nrNodes, nrNodesNew); 
  e4n = getElements4Nodes(n4e);

  uNew = zeros(nrSidesNew, 1);
  val = [u(s4e(:, 1)) - u(s4e(:, 2)) + u(s4e(:, 3)), ... % 1st local nodes
            u(s4e(:, 1)) + u(s4e(:, 2)) - u(s4e(:, 3)), ... % 2nd local nodes
            - u(s4e(:, 1)) + u(s4e(:, 2)) + u(s4e(:, 3))]'; % 3rd local nodes
  val = val(:); 
    % values in the nodes computed by formula above, every three entries are
    % wrt. to one triangle (value in midpoint is (node1-node2)/2, hence 3
    % equations, solve for node1, node 2, and node 3 to obtain above formula)
  valNew = computeP1Extension(n4e, n4eNew, nrElemsNew, val, ...
    n4parentSides4n, e4n); 
    % now possible since values in nodes are knwon
    
    % TODO continue here
  for elem = 1:nrElemsNew % this might be vertorizable somehow, right?
    sides = s4eNew(elem, :);
    val = valNew(3*elem - [2 1 0]); % values in nodes on current element
    uNew(sides) = (val + val([2 3 1]))/2; 
     % middle values of nodes for the current edge for value in midpoint
     % no addition needed since CR continuous in the middle points
     % BUT this overrides already known values for every inner edge with the same value
     % is this fixable?
  end
end

function n4parentSides4n = getParentSide(c4n, c4nNew, s4n, n4sNew, ...
    nrNodes, nrNodesNew)
% return the two nodes of bisected side for new node and twice the input
% for old node
  n4parentSides4n = [[1:nrNodes; 1:nrNodes]'; zeros(nrNodesNew - nrNodes, 2)];
  for node = nrNodes+1:nrNodesNew 
    [row, ~] = find(n4sNew==node);
    n4r = n4sNew(row, :); 
      % n4r - nodes for rows (i.e. sides) in n4s in which node exists (there 
      % are at least two sides since node is a new node)
    list = n4r(n4r<=nrNodes)';
      % list - all nodes in n4r that are oldNodes (at least two (see above))
    if length(list)==2 
      % if there are only two nodes in list then those two must be the side in
      % n4s (without order) of which node is the midpoint of
      n4parentSides4n(node, :) = list;
        % therefore list contains only the nodes of the parent side of node
    else
      com = nchoosek(list, 2); % all possible node combinations of nodes in
        % list i.e. all parent side candidates of node
      mids = computeMid4s(c4n, com);
      I = ismember(mids, c4nNew(node, :), 'rows'); 
        % finds the index of the midpoint in mids that shares its coordinates
        % with node, i.e.  which edge is the parent edge (this is why the
        % prolongation has to be done before the edge nodes are projected
        % outward, else the coordinates wouldn't be correct)
      currParentSide = com(I, :); 
      if size(currParentSide, 1)>1 
        % if there is the more than one candidate for the parent side check
        % which of the nodepairs actually belongs to a side of the coarse
        % triangulation using s4n (this side is the parent side, since its 
        % midpoint is node)
        for k = 1:size(currParentSide, 1)
          if s4n(currParentSide(k, 1), currParentSide(k, 2)) ~= 0
           currParentSide = currParentSide(k, :); 
           break;
          end
        end
      end
      n4parentSides4n(node, :) = currParentSide;
    end
  end
end

function e4n = getElements4Nodes(n4e) 
% find all triangles in the nodepatch of a node
  nrNodes = max(n4e(:));
  e4n = mat2cell(repmat(sparse(nrNodes, nrNodes), nrNodes, 1), ...
    zeros(1, nrNodes) + nrNodes);
    % e4n{j}(k,l) either yields the element with the nodes j, k, l or all zero
    % sparse 1x1 if there is no such element
  for elem = 1:size(n4e, 1)
    P = perms(n4e(elem, :)); % all permutations of the nodes of the element
    for k = 1:size(P, 1)
      e4n{P(k, 1)}(P(k, 2), P(k, 3)) = elem;
    end
  end
end

function valNew = computeP1Extension(n4e,n4eNew,nrElemsNew,vOld,...
                                   n4parentSides4n,e4n)
% set the local values in the new nodes as linear interpolation of the values
% of the vertices of the new nodes parent edge on each triangle (standard nodal
% interpolation of the discrete solution of the old mesh to the new mesh of
% discontinuous P1 functions)
  valNew = zeros(3*nrElemsNew, 1);
  for elem = 1:nrElemsNew
    nodes = n4eNew(elem, :);
    currParentSides = n4parentSides4n(nodes, :); 
      % the nodes of the three parent sides of the nodes of elem
    currParentElem = unique(currParentSides(:)); 
      % nodes of current parent element (without usual order)
    oldElemNumber = e4n{currParentElem(1)}...
      (currParentElem(2), currParentElem(3)); 
      % number of the parent element using e4n
    helper = zeros(3, 2);
    for k = 1:3
      I = currParentSides==n4e(oldElemNumber, k); 
        % get Indices of k-th node of the parent element in currParentSides
      helper(I) = vOld(3*(oldElemNumber - 1) + k); 
        % formula yields the value in the local k-th node of the parent element
        % which is written in helper at the entries specified by I
    end
    valNew(3*elem - [2 1 0]) = sum(helper, 2)/2; 
      % writes the mean value of the values of the two nodes of the parent edge
      % in the current Parent element in the right entry of valNew
  end
end

