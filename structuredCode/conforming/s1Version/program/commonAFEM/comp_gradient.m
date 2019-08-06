function du = comp_gradient(c4n,n4e,u)
nE = size(n4e,1); 
du = zeros(nE,2);
for j = 1:nE
    X_T = [ones(1,3);c4n(n4e(j,:),:)'];
    grads_T = X_T\[zeros(1,2);eye(2)];
      % Since u is S^1 function du=d(u(P1)phi1+u(P2)phi2+u(P3)phi3)
      % = sum(j=1,2,3) ( u(Pj)d(phi_j)) )
    du(j,:) = u(n4e(j,:))'*grads_T;          
end

