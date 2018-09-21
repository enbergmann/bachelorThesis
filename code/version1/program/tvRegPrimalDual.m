function  [u,corr,corr_vec,energy_vec,errorExactVec] = ...
    tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,u,Lambda,f,alpha,epsStop,uExact) 
    % Input red: Number of starting triangulations
    
    %% Initialisation
  
    tau = 1/2;

    nrElems = size(n4e,1);
    area4e = computeArea4e(c4n,n4e);
    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));

    dof = computeDof(n4e,nrSides,n4sDb,n4sNb);

    [STIMANC,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);
    A = STIMANC/tau+alpha*MAMANC; 

    [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);

    du = computeGradientNC(c4n,n4e,u);

    v = zeros(nrSides,1);    

    corr = epsStop+1; 
    corr_vec = [];
    energy_vec = [];
    errorExactVec = [];

    E = 1;
    while corr > epsStop
        dv = computeGradientNC(c4n,n4e,v);
        M = Lambda + tau*(du + tau*dv);
        Lambda = bsxfun(@rdivide,M,max(1,sqrt(sum(M.^2,2))));
        
        [b,temp] = computeRHS(c4n,n4e,s4e,nrSides,area4e,du,tau,Lambda,nrElems,temp1,temp2,temp3);     

        %% Solve System
        uNew = zeros(nrSides,1);
        uNew(dof) = A(dof,dof)\b(dof);
        v=(uNew-u)/tau;        
        u = uNew;

        %% Check Termination
        du = computeGradientNC(c4n,n4e,u);
        ENew = computeEnergy(area4e,u,du,alpha,temp,MAMANC);

        corr = abs(ENew-E);
        fprintf('corr/epsStop: %e / %e \n',corr,epsStop);
        format long;
        fprintf('E = %f, E_exact = %f\n', E, -2.05802391003896);
        format short;
        fprintf('===============================================\n');
        E = ENew;
        energy_vec(end+1) = E;
        errorExactVec(end+1) = sum(error4eCRL2(c4n,n4e,uExact,u)); 
        corr_vec(end+1) = corr;

%     plotCR(c4n,n4e,uNew);
%     clf('reset');
    %     fprintf('\n')
    end
end
