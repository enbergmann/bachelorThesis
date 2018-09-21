function [STIMANC,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems)
  %% Create the stiffness matrix, mass matrix and right-hand side b    
  STIMANClocal = zeros(3,3,nrElems);
  MAMANClocal = zeros(3,3,nrElems);
  for elem = 1 : nrElems
      nodes = n4e(elem,:);   % nodes of this element
      coords = c4n(nodes,:); % coordinates for the nodes
      area = area4e(elem);   % area of this element
      gradsNC = [1 1 1; coords']\[0 0; -2 0; 0 -2]; % gradients for CR basis
      gradsNC = gradsNC([3 1 2],:); % reorder to fit DoF numbering
      STIMANClocal(:,:,elem) = area * (gradsNC * gradsNC'); % local stiffness matrix
      MAMANClocal(:,:,elem) = area * eye(3)/3; % local mass matrix
  end
  
  % assembly of the global stiffness matrix and global mass matrix
  s4eT = s4e';
  I = [s4eT;s4eT;s4eT];
  J = [s4eT(:),s4eT(:),s4eT(:)]';
  STIMANC = sparse(I(:),J(:),STIMANClocal(:));
  MAMANC = sparse(I(:),J(:),MAMANClocal(:));
end
