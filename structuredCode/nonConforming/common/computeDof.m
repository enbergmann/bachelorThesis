function dof = computeDof(n4e,nrSides,n4sDb,n4sNb)
  s4n = computeS4n(n4e);
  
  % Dirichlet boundary sides
  DbSides = zeros(1,size(n4sDb,1));
  for i = 1:size(n4sDb,1)
      DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
  end
      
  % Neumann boundary sides
  NbSides = zeros(1,size(n4sNb,1));
  for i = 1:size(n4sNb,1)
      NbSides(i) = s4n(n4sNb(i,1),n4sNb(i,2));
  end
    
  % degrees of freedom: one per non-Dirichlet side
  dof = setdiff(1:nrSides,DbSides);
end
