function u_fine = computeP1prolongation(c4n,n4e,c4n_fine,u)
  % extend a P1 function to a refined mesh
  n4s = computeN4s(n4e); 
  mid4s = computeMid4s(c4n,n4s);
  c4newNodes = c4n_fine(size(c4n,1)+1:end,:);
  [~,I] = ismember(c4newNodes,mid4s,'rows');
  u_fine = [u; sum(u(n4s(I,:)),2)/2];
end
