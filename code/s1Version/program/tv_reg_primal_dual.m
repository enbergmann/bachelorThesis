function [corrVec,energyVec,c4n,u]=...
  tv_reg_primal_dual(red,terminate,alpha,delta)

  h = 2^(-red); 
  tau = h^(1/2)/10; 
  
  [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red); 
  
  [s,m] = fe_matrices(c4n,n4e);
  ms = mixed_matrix(c4n,n4e);
  A = m+h*s; 
  
  nC = size(c4n,1); 
  nE = size(n4e,1); 
  
  f = @(x)g(x,alpha,delta);  
  uExact = @(x)gUexact(x,alpha,delta);

  % u = zeros(nC,1); 
  u = f(c4n);

  u_tilde = u; 
  p = zeros(nE,2); 

  DbNodes = unique(n4sDb);    
  dof = setdiff(1:nC,DbNodes);

  corr = 1; % ||d_t u_h^k||_{h,1/2} 
  while corr > terminate
      % step=step+1;
      du_tilde = comp_gradient(c4n,n4e,u_tilde);
      p_tmp = p+tau*du_tilde;
      p = p_tmp./max(1,(sqrt(sum(p_tmp.^2,2))*ones(1,2))); 
      P = reshape(p',2*nE,1);
      
      u_new = zeros(nC,1);
      LHS = (A+tau*alpha*m);
      RHS = (A*u-tau*ms*P+tau*alpha*m*f(c4n));
      u_new(dof) = LHS(dof,dof)\RHS(dof);
      dt_u = (u-u_new)/tau; 
      u_tilde = 2*u_new-u; 
      
      % Bartels termination criterion
      corr = sqrt(dt_u'*A*dt_u);

      fprintf('corr/epsStop: %e / %e\n',corr,terminate);

      u = u_new;
      show_p1(c4n,n4e,n4sDb,n4sNb,u);
  end
end
