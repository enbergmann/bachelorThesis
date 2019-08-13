function du = gradientCR(c4n,n4e,u)
% Computes the piecewise gradient of the Crouzeix-Raviart function u (vanishing
% in the midpoints of boundary edges) with respect to the triangulation given
% by [c4n, n4e].
%
% gradientCR.m
% input:  c4n - coordinates for nodes
%         n4e - nodes for elements
%         u   - 'function handle' of the function whose piecewise gradient is
%               to be computed
%
% output: du  - (number of elements x 2)-dimensional array where the j-th row
%               contains the gradient of u on the j-th triangle 

  nrElems = size(n4e,1); 
  s4e = computeS4e(n4e);
  du = zeros(nrElems,2);
  for elem = 1:nrElems
      gradsT = [ones(1,3); c4n(n4e(elem,:),:)']\[zeros(1,2);-2*eye(2)];
      gradsT = grads_T([3 1 2],:);
      du(elem,:) = u(s4e(elem,:))'*gradsT;          
  end
end
