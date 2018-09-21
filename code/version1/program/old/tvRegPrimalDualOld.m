function  [u,corr,corr_vec,energy_vec] = tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,u,Lambda,f,alpha,epsStop) 
    % Input red: Number of starting triangulations
    
    %% Initialisation
  
    tau = 1/2;

    nrElems = size(n4e,1);
    area4e = computeArea4e(c4n,n4e);
    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));

    dof = computeDof(n4e,nrSides,n4sDb,n4sNb);

          %% Create the stiffness matrix, mass matrix and right-hand side b    
      STIMANClocal = zeros(3,3,nrElems);
      MAMANClocal = zeros(3,3,nrElems);
      for elem = 1 : nrElems
          nodes = n4e(elem,:);   % nodes of this element
          coords = c4n(nodes,:); % coordinates for the nodes
          area = area4e(elem);   % area of this element
          gradsNC = [1 1 1; coords']\[0 0; -2 0; 0 -2]; % gradients for CR basis
          gradsNC = gradsNC([3 1 2],:); % reorder to fit DoF numbering
          STIMANClocal(:,:,elem) = area * (gradsNC * gradsNC'); % local stiffness matrix
          MAMANClocal(:,:,elem) = area * eye(3)/3; % local mass matrix
      end

      % assembly of the global stiffness matrix and global mass matrix
      s4eT = s4e';
      I = [s4eT;s4eT;s4eT];
      J = [s4eT(:),s4eT(:),s4eT(:)]';
      STIMANC = sparse(I(:),J(:),STIMANClocal(:));
      MAMANC = sparse(I(:),J(:),MAMANClocal(:));
        A = STIMANC/tau+alpha*MAMANC; 
    
    degree = 20;
      phi1 = @(x)( x(:,1) );
      phi2 = @(x)( x(:,2) );
      phi3 = @(x)( 1 - x(:,1) - x(:,2) );

      psi1 = @(x)( 1-2*phi3(x) );
      psi2 = @(x)( 1-2*phi1(x) );
      psi3 = @(x)( 1-2*phi2(x) );

      INT1 = @(n4p,Gpts4p,Gpts4ref) ...
          psi1(Gpts4ref).*f(Gpts4p);
      INT2 = @(n4p,Gpts4p,Gpts4ref) ...
          psi2(Gpts4ref).*f(Gpts4p);        
      INT3 = @(n4p,Gpts4p,Gpts4ref) ...
          psi3(Gpts4ref).*f(Gpts4p);

      temp1 = integrate(INT1, c4n, n4e, degree, area4e);
      temp2 = integrate(INT2, c4n, n4e, degree, area4e);
      temp3 = integrate(INT3, c4n, n4e, degree, area4e);

    du = computeGradientNC(c4n,n4e,u);

    v = zeros(nrSides,1);    

    corr = epsStop+1; 
    corr_vec = [];
    energy_vec = [];

    E = 1;
    while corr > epsStop  
          dv = computeGradientNC(c4n,n4e,v);
          M = Lambda + tau*(du + tau*dv);
          Lambda = bsxfun(@rdivide,M,max(1,sqrt(sum(M.^2,2))));
          %% Create the right-hand side b  (midpoint rule might be seperated from integrate, integrate wo loop)
          b = zeros(nrSides,1);
          for elem = 1 : nrElems
              nodes = n4e(elem,:);   % nodes of this element
              sides = s4e(elem,:);   % sides of this element
              coords = c4n(nodes,:); % coordinates for the nodes
              area = area4e(elem);   % area of this element
              gradsNC = [1 1 1; coords']\[0 0; -2 0; 0 -2]; % gradients for CR basis
              gradsNC = gradsNC([3 1 2],:); % reorder to fit DoF numbering

              bLocal = (du(elem,:)/tau - Lambda(elem,:)) * gradsNC';  
              b(sides) = b(sides) + [temp1(elem),temp2(elem),temp3(elem)]' +area*bLocal'; % right-hand side
              % midpoint rule
        %       mid = mid4e(elem,:);     % midpoint of this element
        %       b(sides) = b(sides) + area*f(mid)*ones(3,1)/3 +area*bLocal'; % right-hand side
          end

        %% Solve System
        uNew = zeros(nrSides,1);
        uNew(dof) = A(dof,dof)\b(dof);
        v=(uNew-u)/tau;        
        u = uNew;

        %% Check Termination
        du = computeGradientNC(c4n,n4e,u);
        ENew = computeEnergy(nrElems,s4e,area4e,u,du,alpha,temp1,temp2,temp3);

%        corr = norm(v)/norm(u);
        corr = abs(ENew-E);
        fprintf('corr/epsStop: %f / %f \n',corr,epsStop)
        fprintf('E --> Enew: %f --> %f \n', E, ENew)
        E = ENew;
        energy_vec(end+1) = E;
        corr_vec(end+1) = corr;

     plotCR(c4n,n4e,uNew);
     clf('reset');
    %     fprintf('\n')
    end
end

%% Functions

function dof = computeDof(n4e,nrSides,n4sDb,n4sNb)
  s4n = computeS4n(n4e);
  
  % Dirichlet boundary sides
  DbSides = zeros(1,size(n4sDb,1));
  for i = 1:size(n4sDb,1)
      DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
  end
      
  % Neumann boundary sides
  NbSides = zeros(1,size(n4sNb,1));
  for i = 1:size(n4sNb,1)
      NbSides(i) = s4n(n4sNb(i,1),n4sNb(i,2));
  end
    
  % degrees of freedom: one per non-Dirichlet side
  dof = setdiff(1:nrSides,DbSides);
end


function ENew = computeEnergy(nrElems,s4e,area4e,u,du,alpha,temp1,temp2,temp3)
  ENew = 0;
    for elem = 1 : nrElems
%       nodes = n4e(elem,:);   % nodes of this element
      sides = s4e(elem,:);   % sides of this element
%       coords = c4n(nodes,:); % coordinates for the nodes
      area = area4e(elem);   % area of this element
%       gradsNC = [1 1 1; coords']\[0 0; -2 0; 0 -2]; % gradients for CR basis
%       gradsNC = gradsNC([3 1 2],:); % reorder to fit DoF numbering
      coeff = u(sides); % coefficients on this element
      
      q = sum(coeff.^2); % sum of quadratic terms
      s = sum(coeff);
      m = coeff([2 3 1])'*coeff; % sum of mixed terms

%       ENew = ENew + area*( ...
%         + alpha/6*( -s^2 + 2*q + 2*m ) ... % L^2 norm of u
%         + norm(du(elem)) ) ... % L^2 norm of gradient(NC) of u
%         - coeff'*[temp1(elem),temp2(elem),temp3(elem)]'; % integral over u*f
    
      ENew = ENew + area*( ...
        + alpha/6*( -q+2*m+2*(q+m) ) ... % L^2 norm of u
        + norm(du(elem)) ) ... % L^2 norm of gradient(NC) of u
        - coeff'*[temp1(elem),temp2(elem),temp3(elem)]'; % integral over u*f
    end
end



