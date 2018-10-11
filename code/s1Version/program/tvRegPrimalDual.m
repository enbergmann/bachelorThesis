function [corrVec,energyVec,u]=...
  tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,terminate,alpha,f,u)
  
  [s,m] = fe_matrices(c4n,n4e);
  ms = mixed_matrix(c4n,n4e);
  A = m+h*s; 
  
  nC = size(c4n,1); 
  nE = size(n4e,1); 
  
  u_tilde = u; 
  p = zeros(nE,2); 

  DbNodes = unique(n4sDb);    
  dof = setdiff(1:nC,DbNodes);

  corr = 1; % ||d_t u_h^k||_{h,1/2} 

  area4e = computeArea4e(c4n,n4e);
  s4e = computeS4e(n4e);
  nrSides = max(max(s4e));
  [~,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nE);
  n4s = computeN4s(n4e);
  [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);
  temp = zeros(nrSides,1);
  for elem = 1:1:nE
    temp(s4e(elem,:)) = temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
  end
  corrVec=[];
  energyVec=[];
  E = 1;

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
      
      corr = sqrt(dt_u'*A*dt_u); % Bartels termination criterion
      
      uNC = convertS1toCR(n4s,u);
      duNC = computeGradientNC(c4n,n4e,uNC);
      ENew = computeEnergy(area4e,uNC,duNC,alpha,temp,MAMANC);
      E = ENew;

      fprintf('corr/epsStop: %e / %e\n',corr,terminate);
      format long;
      fprintf('E = %f, E_exact = %f\n', E, -2.05802391003896);
      format short;
      fprintf('============================== \n');

      corrVec=cat(2,corrVec,corr);
      
      energyVec=cat(2,energyVec,ENew);

      u = u_new;
      show_p1(c4n,n4e,u);
  end
end
