function [stiMaNC, maMaNC] = computeFeMatricesCR(currData)
  %% Create the stiffness matrix, mass matrix and right-hand side b    
  %c4n,n4e,s4e,area4e,nrElems
  
  %TODO
% Computes the gradients of the Crouzeix-Raviart basis functions on every
% triangle with respect to the triangulation given by [c4n, n4e].
%
% computeFeMatricesCR.m
% input:  currData - struct with fields:
%                        nrElems: number of elements
%                         area4e: areas for elements
%                      gradsCR4e: gradients of CR basis functions of elements
%                            s4e: sides for elements
%
% output: stiMaNC  - '(3 x 2 x nrElems)-dimensional double array' where the 
%                    j-th row of the k-th matrix contains the gradient of the
%                    CR-basis function w.r.t. to j-th edge of the k-th element
%         maMaNC   -
  
  % extract necessary data
  nrElems = currData.nrElems;
  area4e = currData.area4e;
  gradsCR4e = currData.gradsCR4e;
  s4e = currData.s4e;

  % compute local stiffness and mass matrices
  stiMaNClocal = zeros(3, 3, nrElems);
  maMaNClocal = zeros(3, 3, nrElems);
  for elem = 1:nrElems
    stiMaNClocal(:, :, elem) = ...
      area4e(elem)*(gradsCR4e(:, :, elem)*gradsCR4e(:, :, elem)'); 
      % local stiffness matrix
    maMaNClocal(:, :, elem) = area4e(elem)*eye(3)/3; % local mass matrix
  end
  
  % assembly of the global stiffness matrix and global mass matrix
  s4eT = s4e';
  I = [s4eT; s4eT; s4eT];
  J = [s4eT(:), s4eT(:), s4eT(:)]';
  stiMaNC = sparse(I(:), J(:), stiMaNClocal(:));
  maMaNC = sparse(I(:), J(:), maMaNClocal(:));
end
