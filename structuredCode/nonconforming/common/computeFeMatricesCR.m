function [stiMaCR, maMaCR] = computeFeMatricesCR(currData)
%% DOC
% Assembles the Crouzeix-Raviart stiffness and mass matrix with respect to the
% triangulation given by [c4n, n4e].
% 
% computeFeMatricesCR.m
% input:  currData - 'struct' with fields:
%                        nrElems: number of elements
%                         area4e: areas for elements
%                      gradsCR4e: '(3 x 2 x nrElems)-dimensional double
%                                 array' where the j-th row of the k-th matrix
%                                 contains the gradient of the CR-basis
%                                 function w.r.t. to j-th edge of the k-th
%                                 element
%                            s4e: sides for elements
%
% output: stiMaCR  - '(nrSides x nrSides)-dimensional sparse double array'
%                    where the k-th entry in the j-th row is the L2-scalar
%                    product of the piecewise gradient of the CR-basis function
%                    w.r.t. k-th edge with the piecewise gradient of the
%                    CR-basis function w.r.t. j-th edge
%         maMaCR   - '(nrSides x nrSides)-dimensional sparse double array'
%                    where the k-th entry in the j-th row is the L2-scalar
%                    product of the CR-basis function w.r.t. k-th edge with the
%                    CR-basis function w.r.t. j-th edge
  
%% INIT
  % extract necessary information from currData
  nrElems = currData.nrElems;
  area4e = currData.area4e;
  gradsCR4e = currData.gradsCR4e;
  s4e = currData.s4e;

  % initialize local stiffness and mass matrices
  stiMaCRlocal = zeros(3, 3, nrElems);
  maMaCRlocal = zeros(3, 3, nrElems);

%% MAIN
  % compute local stiffness and mass matrices
  for elem = 1:nrElems
    stiMaCRlocal(:, :, elem) = ...
      area4e(elem)*(gradsCR4e(:, :, elem)*gradsCR4e(:, :, elem)'); 
      % local stiffness matrix
    maMaCRlocal(:, :, elem) = area4e(elem)*eye(3)/3; 
      % local mass matrix
  end
  
  % assembly of the global stiffness matrix and global mass matrix
  s4eT = s4e';
  I = [s4eT; s4eT; s4eT];
  J = [s4eT(:), s4eT(:), s4eT(:)]';
  stiMaCR = sparse(I(:), J(:), stiMaCRlocal(:));
  maMaCR = sparse(I(:), J(:), maMaCRlocal(:));
end
