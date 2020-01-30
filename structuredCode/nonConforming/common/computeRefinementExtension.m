% TODO
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

% function [uNew] = computeCRprolongation(c4n, n4e, c4nNew, n4eNew, u)
function [uNew] = computeRefinementExtension(c4n, n4e, c4nNew, n4eNew, u)

  %% initialization
  n4sNew = computeN4s(n4eNew);
  s4n = computeS4n(n4e);
  s4e = computeS4e(n4e); 
  s4eNew = computeS4e(n4eNew);

  nrNodes = size(c4n, 1); 
  nrNodesNew = size(c4nNew, 1);
  nrElemsNew = size(n4eNew, 1);
  nrSidesNew = size(n4sNew, 1);

  %TODO continue here use mlint
  n4parentSides4n = getParentSide(1:nrNodesNew,c4n,c4nNew,s4n,n4sNew,nrNodes); 
  %nodes for parent side for every node
  e4nOld = getElements4Nodes(n4e);


      uNew = zeros(nrSidesNew,1);
      valOld = [u(s4e(:,1))-u(s4e(:,2))+u(s4e(:,3)),...  
                u(s4e(:,1))+u(s4e(:,2))-u(s4e(:,3)),...
                - u(s4e(:,1))+u(s4e(:,2))+u(s4e(:,3))]';
      valOld = valOld(:); % values in the nodes computed by formula above, every three entries are wrt to one triangle
      % value in middlepoint is (node1-node2)/2, hence 3 equations, solve for node1, 2, 3, yields above result
      valNew = computeP1Extension(n4e,n4eNew,nrElemsNew,valOld,... 
                                  n4parentSides4n,e4nOld); % now possible since values in nodes are knwon
      for elem = 1:nrElemsNew % this might be vertorizable somehow, right?
        sides = s4eNew(elem,:);
        val = valNew(3*elem-[2 1 0]); % values in nodes on current element
        uNew(sides) = (val+val([2 3 1]))/2; % middle values of nodes for the current edge for value in middlepoint
         % no addition needed since CR continuous in the middle points
         % BUT this overrides already known values for every inner edge with the same value
         % is this fixable?
      end
end

% set the values in the new nodes as linear interpolation of the values of
% the vertices of the new nodes parent edge on each triangle
% (standard nodal interpolation of the discrete solution of the old mesh to
% the new mesh of discontinuous P1 functions)
function vNew = computeP1Extension(n4e,n4eNew,nrElemsNew,vOld,...
                                   n4parentSides4n,e4nOld)
  vNew = zeros(3*nrElemsNew,1);
  for elem = 1:nrElemsNew
    nodes = n4eNew(elem,:);
    currParentSide = n4parentSides4n(nodes,:); % parent sides of all nodes in elem
    currParentElem = unique(currParentSide(:)); % nodes of current parent element (without usual order)
    loc = e4nOld{currParentElem(1)};
    oldElemNumber = loc(currParentElem(2),currParentElem(3)); % this and the line above could just be one line, (for real) e4nOld{currParentElem(1)}(currParentElem(2),currParentElem(3))
    helper = zeros(size(currParentSide));
    for k=1:3
      I = find(currParentSide==n4e(oldElemNumber,k)); 
      %get Indices (not subscripts) of k-th node of the parent element in the parent sides
      helper(I) = vOld(3*(oldElemNumber-1)+k); % formula yields the k-th node in the parent element, hence helper gets those values in the right place given by I
    end
    vNew(3*elem-[2 1 0]) = sum(helper,2)/2; %middles the values of the two nodes of the parent edge in the current Parent element and writes it at the right place in vNew
  end
end

% return the two nodes of bisected side for new node and twice the input
% for old node
function n4parentSides4n = getParentSide(nodeList,c4n,c4nNew,s4n,...
                                    n4sNew,nrNodes)
% give nrNodesNew instead of nodeList = 1:nrNodesNew

  n4parentSides4n = zeros(length(nodeList),2);
  %nodeList is always 1:nrNodesNew, so pretty unnessacary maybe
  for j = 1:length(nodeList) % for j = 1:nrNodesNew
    node = nodeList(j); % see above, node = j always, so 
    if node <= nrNodes
      n4parentSides4n(j,:) = [node node];
    else
      [row, ~] = find(n4sNew==node);
      n4r = n4sNew(row,:); list = n4r(n4r<=nrNodes)'; %n4r all rows in n4s/ all sides where node is, list all nodes in n4e that are oldNodes
      currParentSide = list(1:2); %all parent sides for node
      if length(list)>2 % if theres only 2 than those two must be the side in n4s (without order), if theres 3 than middle of an old triangle (?)
        com = nchoosek(list,2); mids = computeMid4s(c4n,com);
        I = ismember(mids,c4nNew(node,:),'rows'); % this is why the prolongation
          % has to be done before the edge nodes are projected outward. Which midpoint is the new node
        currParentSide = com(I,:); %parentside is the one node is the midpoit of
        if size(currParentSide,1)>1
          for k=1:size(currParentSide,1)
            if s4n(currParentSide(k,1),currParentSide(k,2)) ~= 0
             currParentSide = currParentSide(k,:); break;
            end
          end
        end
      end
      n4parentSides4n(j,:) = currParentSide;
    end
  end
end

% find all triangles in the nodepatch of a node
function e4n = getElements4Nodes(n4e) % e4n{j}(k,l) either yields the element with the nodes j,k,l or All zero sparse 1x1 if there is no such element
  nrNodes = max(n4e(:));
  e4n = cell(nrNodes,1); for j=1:nrNodes; e4n{j}=sparse(nrNodes,nrNodes); end
  for elem = 1:size(n4e,1)
    nodes = n4e(elem,:); % nodes of element
    P = perms(nodes); % all permutations of nodes
    for k=1:size(P,1)
      e4n{P(k,1)}(P(k,2),P(k,3)) = elem;
    end
  end
end
