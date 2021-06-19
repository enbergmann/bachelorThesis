function vNew = prolongationJ1(c4n, n4e, n4sB, c4nNew, n4eNew, v)
%% DOC
% Computes the prolongation of a CR function v with respect to a mesh defined
% by [c4n, n4e] to a one-level refinement of this mesh defined by [c4nNew,
% n4eNew] using the enriching operator defined by
% J_1 v (z) := |\Tcal(z)|^{-1} \sum_{T\in\Tcal(z)}v|_T(z) 
% for all z\in\Ncal(\Omega).
%
% prolongationJ1.m
% input: c4n    - coordinates for nodes of the coarse mesh
%        n4e    - nodes for elements of the coarse mesh
%        n4sB   - nodes for boundary sides of the coarse mesh
%        c4nNew - coordinates for nodes of the refined mesh
%        n4eNew - nodes for elements of the refined mesh
%        v      - '(nrSides of the coarse mesh x 1)-dimensional double array'
%                 where the j-th row contains the coefficient of the CR
%                 function w.r.t. the j-th side of the coarse triangulation
%
% output: vNew   - '(nrSides of the refined mesh x 1)-dimensional double array'
%                  where the j-th row contains the coefficient of the CR
%                  prolongation of v w.r.t. the j-th side of the new
%                  triangulation

%% INIT
  s4n = computeS4n(n4e);
  n4sNew = computeN4s(n4eNew);
  s4e = computeS4e(n4e); 
  s4eNew = computeS4e(n4eNew);

  nrNodes = size(c4n, 1); 
  nrNodesNew = size(c4nNew, 1);
  nrElemsNew = size(n4eNew, 1);
  nrSidesNew = size(n4sNew, 1);

  n4parentSides4n = getParentSide(c4n, c4nNew, s4n, n4sNew, ...
    nrNodes, nrNodesNew); 
  e4n = getElements4Nodes(n4e);

%% MAIN
  % compute enriching operator
  nodeValuesCR4e = computeNodeValuesCR4e(s4e, v);
  innerNodes = setdiff(1:nrNodes, unique(n4sB(:)));

  nodeValuesJ1 = zeros(size(nodeValuesCR4e));
  for j = innerNodes
    indNode = find(n4e == j);
    localVals4node = nodeValuesCR4e(indNode);
    nodeValuesJ1(indNode) = sum(localVals4node)/length(localVals4node);
  end
  % TODO here probably just 
  % vJ1 = computeJ1(n4e, n4sB, v);
  % nodeValuesJ1 = vJ1(n4e);
  % TODO test if this produces the same here at some point
  % and then it actually probably just is
  % vJ1 = computeJ1(n4e, n4sB, v);
  % v = courant2CR(...)
  % 
  % thats not true, computeP1 extension is still inbetween
  
  % compute prolongation using enriching operator
  val = transpose(nodeValuesJ1); 
  val = val(:);
    % every three entries are w.r.t. to one triangle   
  valNew = computeP1Extension(n4e, n4eNew, nrElemsNew, val, ...
  n4parentSides4n, e4n); 
    % now possible since values in nodes are known
    
  vNew = zeros(nrSidesNew, 1);
  for elem = 1:nrElemsNew
    sides = s4eNew(elem, :);
    temp = valNew(3*elem - [2 1 0]); % values in nodes of current element
    vNew(sides) = (temp + temp([2 3 1]))/2; 
     % average the values in nodes of the current edge for the value in
     % midpoint 
     % no addition needed since CR continuous in the midpoints BUT this
     % overrides already known values for every inner edge with the same value,
     % which is fine since the enriching operator is globally continuous
  end

  % TODO here probably just 
  % see above
  % TODO test if this produces the same here at some point with many steps 
  % inbetween
  % test
  % vJ1 = computeJ1(n4e, n4sB, v);
  % v = courant2CR(...)
  % against the complete old code above
  % 
  % thats not true, computeP1 extension is still inbetween

end

%% FURTHER FUNCTIONS
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

function valNew = computeP1Extension(n4e, n4eNew, nrElemsNew, vOld, ...
                                   n4parentSides4n, e4n)
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
