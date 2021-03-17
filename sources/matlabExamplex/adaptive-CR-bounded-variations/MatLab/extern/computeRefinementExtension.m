% function computeRefinementExtension
%
% Compute the prolongation of a solution of the CRFEM to a
% one-level refinement of the mesh.
%
% INPUT         c4nOld          Coordinates of the coarse mesh
%               n4eOld          nodes for elements of the coarse mesh
%               c4nNew          Coordinates of the fine mesh
%               n4eNew          nodes for elements of the fine mesh
%               solOld          solution on the coarse mesh
%
% OUTPUT        solNew          prolongation of the old solution to the new
%                               mesh

function [solNew] = computeRefinementExtension(c4nOld,n4eOld,c4nNew,...
                           n4eNew,solOld)
  %% initialization
  n4sOld = computeN4s(n4eOld); n4sNew = computeN4s(n4eNew);
  s4nOld = computeS4n(n4eOld);
  s4eOld = computeS4e(n4eOld); s4eNew = computeS4e(n4eNew);

  nNodesOld = size(c4nOld,1); nNodesNew = size(c4nNew,1);
  nElemOld = size(n4eOld,1); nElemNew = size(n4eNew,1);
  nSidesOld = size(n4sOld,1); nSidesNew = size(n4sNew,1);

  nodesNew = [nNodesOld:nNodesNew]';
  % NOTE nodes are not sorted any longer
  sidesNew = setdiff(sort(n4sNew,2),sort(n4sOld,2),'rows');

  parentSide = getParentSide(1:nNodesNew,c4nOld,c4nNew,s4nOld,n4sNew,nNodesOld);
  e4nOld = getElements4Nodes(n4eOld);


      solNew = zeros(nSidesNew,1);
      valOld = [solOld(s4eOld(:,1))-solOld(s4eOld(:,2))+solOld(s4eOld(:,3)),...  
                solOld(s4eOld(:,1))+solOld(s4eOld(:,2))-solOld(s4eOld(:,3)),...
                - solOld(s4eOld(:,1))+solOld(s4eOld(:,2))+solOld(s4eOld(:,3))]';
      valOld = valOld(:);
      
      valNew = computeP1Extension(n4eOld,n4eNew,nElemNew,valOld,...
                                  parentSide,e4nOld);
      for elem = 1:nElemNew
        nodes = n4eNew(elem,:);
        sides = s4eNew(elem,:);
        val = valNew(3*elem-[2 1 0]);
        solNew(sides) = (val+val([2 3 1]))/2;
      end

end

% set the values in the new nodes as linear interpolation of the values of
% the vertices of the new nodes parent edge on each triangle
% (standard nodal interpolation of the discrete solution of the old mesh to
% the new mesh of discontinuous P1 functions)
function vNew = computeP1Extension(n4eOld,n4eNew,nElemNew,vOld,...
                                   parentSide,e4nOld)
  vNew = zeros(3*nElemNew,1);
  for elem = 1:nElemNew
    nodes = n4eNew(elem,:);
    currParentSide = parentSide(nodes,:);
    currParentElem = unique(currParentSide(:));
    loc = e4nOld{currParentElem(1)};
    oldElemNumber = loc(currParentElem(2),currParentElem(3));
    helper = zeros(size(currParentSide));
    for k=1:3
      I = find(currParentSide==n4eOld(oldElemNumber,k));
      helper(I) = vOld(3*(oldElemNumber-1)+k);
    end
    vNew(3*elem-[2 1 0]) = sum(helper,2)/2;
  end
end

% return the two nodes of bisected side for new node and twice the input
% for old node
function parentSide = getParentSide(nodeList,c4nOld,c4nNew,s4nOld,...
                                    n4sNew,nNodesOld)
  parentSide = zeros(length(nodeList),2);
  for j = 1:length(nodeList)
    node = nodeList(j);
    if node <= nNodesOld
      parentSide(j,:) = [node node];
    else
      [row,col] = find(n4sNew==node);
      n4r = n4sNew(row,:); list = n4r(n4r<=nNodesOld)';
      currParentSide = list(1:2);
      if length(list)>2
        com = nchoosek(list,2); mids = computeMid4s(c4nOld,com);
        I = ismember(mids,c4nNew(node,:),'rows');
        currParentSide = com(I,:);
        if size(currParentSide,1)>1
          for k=1:size(currParentSide,1)
            if s4nOld(currParentSide(k,1),currParentSide(k,2))~=0;
             currParentSide = currParentSide(k,:); break;
            end
          end
        end
      end
      parentSide(j,:) = currParentSide;
    end
  end
end

% find all triangles in the nodepatch of a node
function e4n = getElements4Nodes(n4e)
  nNodes = max(n4e(:));
  e4n = cell(nNodes,1); for j=1:nNodes; e4n{j}=sparse(nNodes,nNodes); end
  for elem = 1:size(n4e,1)
    nodes = n4e(elem,:);
    P = perms(nodes);
    for k=1:size(P,1)
      e4n{P(k,1)}(P(k,2),P(k,3)) = elem;
    end
  end
end
