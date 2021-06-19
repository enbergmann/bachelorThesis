%TODO
%can probably use the function below, just that parentsside is
%every side or sth

% instead of nodeValuesJ1 we just in general have for the Courant fct v
% nodeValues = v(n4e)
% and can then use the rest of the code as usual
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
