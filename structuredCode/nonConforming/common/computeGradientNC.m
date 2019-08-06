function du = compGradientNC(c4n,n4e,u)
  nrElems = size(n4e,1); 
  s4e = computeS4e(n4e);
  du = zeros(nrElems,2);
  for elem = 1:nrElems
      grads_T = [ones(1,3);c4n(n4e(elem,:),:)']\[zeros(1,2);-2*eye(2)];
      grads_T = grads_T([3 1 2],:);
      du(elem,:) = u(s4e(elem,:))'*grads_T;          
  end
end
