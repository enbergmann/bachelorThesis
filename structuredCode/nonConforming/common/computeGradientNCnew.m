function du = computeGradientNCnew(c4n,n4e,u,grads4e)
  nrElems = size(n4e,1); 
  s4e = computeS4e(n4e);
  du = zeros(nrElems,2);
  for elem = 1:nrElems
      du(elem,:) = u(s4e(elem,:))'*grads4e(:,:,elem);          
  end
end
