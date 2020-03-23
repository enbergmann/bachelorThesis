function dof = computeDofCR(currData)
%% DOC
% Computes the degrees of freedom for the CR_0^1 FEM of the triangulation given
% by [n4e, n4sDb, n4sNb].
%
% computeDofCR.m
% input: currData - 'struct' with fields:
%                       n4sDb: nodes for Dirichlet boundary sides
%                         s4n: sides for nodes
%                     nrSides: number of sides
%
% output: dof - '(1 x nrDof)-dimensional double array' where the j-th row
%               contains the number of the j-th degree of freedom
  
%% INIT
  % extract necessary information from currData
  n4sDb = currData.n4sDb;
  s4n = currData.s4n;
  nrSides = currData.nrSides;

%% MAIN
  % compute Dirichlet boundary sides
  DbSides = zeros(1, size(n4sDb, 1));
  for j = 1:size(n4sDb, 1), DbSides(j) = s4n(n4sDb(j, 1), n4sDb(j, 2)); end
    % get the side number of the j-th Dirichlet boundary side via its nodes
    
  % compute degrees of freedom: one per non-Dirichlet side
  dof = setdiff(1:nrSides, DbSides);
end
