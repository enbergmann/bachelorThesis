function gradsCR4e = computeGradsCR4e(currData)
%% DOC
% Computes the gradients of the Crouzeix-Raviart basis functions on every
% triangle with respect to the triangulation given by [c4n, n4e].
%
% computeGradsCR4e.m
% input: currData - 'struct' with fields:
%                         c4n: coordinates for nodes
%                         n4e: nodes for elements
%                     nrElems: number of elements
%
% output: gradsCR4e - '(3 x 2 x nrElems)-dimensional double array' where the 
%                     j-th row of the k-th matrix contains the gradient of the
%                     CR-basis function w.r.t. the j-th edge of the k-th element
  
%% INIT
  % extract necessary data
  nrElems = currData.nrElems;
  c4n = currData.c4n;
  n4e = currData.n4e;

  % initialize gradients
  gradsCR4e = zeros(3, 2, nrElems);

%% MAIN
  % compute gradsCR4e
  for elem = 1:nrElems
    gradsT = [ones(1, 3); c4n(n4e(elem, :), :)']\[zeros(1, 2); -2*eye(2)];
    gradsT = gradsT([3 1 2], :);
    gradsCR4e(:, :, elem) = gradsT;
  end
end
