function gradCRv = gradientCR(currData, v)
% Computes the piecewise gradient of the Crouzeix-Raviart function v (vanishing
% in the midpoints of boundary edges) with respect to the triangulation given
% by [c4n, n4e].
%
% gradientCR.m
% input:  currData - struct with fields:
%                  nrElems: number of elements
%                      s4e: sides for elements
%                gradsCR4e: gradients of side based Crouzeix-Raviart basis functions
%                           for all elements
%         v        - 'function_handle' of the function whose piecewise gradient 
%                     is to be computed
%
% output: gradCRv  - '(number of elements x 2)-dimensional double array' where 
%                     the j-th row contains the gradient of v on the j-th
%                     triangle 

  % extract necessary data
  nrElems = currData.nrElems;
  s4e = currData.s4e;
  gradsCR4e = currData.gradsCR4e;

  % compute uCR
  gradCRv = zeros(nrElems,2);
  for elem = 1:nrElems
    gradCRv(elem,:) = v(s4e(elem,:))'*gradsCR4e(:,:,elem);          
  end
end
