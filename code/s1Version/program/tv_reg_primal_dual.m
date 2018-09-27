function [c4n,n4e,u]=tv_reg_primal_dual(red,terminate)

  alpha = 1;
  delta = 1;
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

  u = zeros(nC,1); 

  % u = f(c4n);
  u_tilde = u; 
  p = zeros(nE,2); 

  % nrElems = size(n4e,1);
  % area4e = computeArea4e(c4n,n4e);
  % s4e = computeS4e(n4e);
  % nrSides = max(max(s4e));
  % [~,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);
  % n4s = computeN4s(n4e);

  DbNodes = unique(n4sDb);    
  dof = setdiff(1:nC,DbNodes);
  % x = zeros(nC,1); % get the Dirichlet values
  % DbCoords = c4n(DbNodes,:); % coordinates of Dirichlet nodes
  % x(DbNodes) = u4Db(DbCoords);
  % b = b - A * x;  % substract inhomogenous boundary


  % [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);
  % temp = zeros(nrSides,1);
  % for elem = 1:1:nrElems
  %   temp(s4e(elem,:)) = temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
  % end

  corr = 1; % ||d_t u_h^k||_{h,1/2} 
  terminate = 1e-2;
  % corr_vec=[];
  % energy_vec=[];
  % errorExactVec=[];
  % step=0;
  % E = 1;
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
      % u_new = LHS\RHS;
      dt_u = (u-u_new)/tau; 
      u_tilde = 2*u_new-u; 
      
      corr = sqrt(dt_u'*A*dt_u);

      % uNC = convertS1toCR(n4s,u);
      % duNC = computeGradientNC(c4n,n4e,uNC);
      % ENew = computeEnergy(area4e,uNC,duNC,alpha,temp,MAMANC);
      % corr = abs(ENew-E);
      % E = ENew;

      u = u_new;
      show_p1(c4n,n4e,n4sDb,n4sNb,u);
  end
end
