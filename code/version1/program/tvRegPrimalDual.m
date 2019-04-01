function  [u,corrVec,energyVec,nrDof] = ...
  tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,epsStop,alpha,f,u,Lambda, ...
                  saveScreenshots) 

  %saveScreenshots: integer x bigger (save results every x iterations)
  %                 or equal to zero (don't save screenshots)

  firstScreenshot = datestr(now,'yy_mm_dd_HH_MM_SS');

  if nargin<13
    saveScreenshots = 0;
  end

  nrElems = size(n4e,1);
  area4e = computeArea4e(c4n,n4e);
  s4e = computeS4e(n4e);
  nrSides = max(max(s4e));

  dof = computeDof(n4e,nrSides,n4sDb,n4sNb);
  nrDof = length(dof);

  [STIMANC,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);
  A = STIMANC/tau+alpha*MAMANC; %TODO here could be an h in front of STIMANC
  C = MAMANC + h*STIMANC;

  [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);

  du = computeGradientNC(c4n,n4e,u);

  v = zeros(nrSides,1);    

  corr = epsStop+1; 
  corrVec = [];
  energyVec = [];
  E = 1;



  nrElems = size(n4e,1); 
  s4e = computeS4e(n4e);
  grads4e = zeros(3,2,nrElems);
  for elem = 1:nrElems
      grads_T = [ones(1,3);c4n(n4e(elem,:),:)']\[zeros(1,2);-2*eye(2)];
      grads_T = grads_T([3 1 2],:);
      grads4e(:,:,elem) = grads_T;
  end



  
  while corr > epsStop
      dv = computeGradientNCnew(c4n,n4e,v,grads4e);
      M = Lambda + tau*(du + tau*dv);
      Lambda = bsxfun(@rdivide,M,max(1,sqrt(sum(M.^2,2))));
      
      [b,temp] = computeRHS(c4n,n4e,s4e,nrSides,area4e, ...
        du,tau,Lambda,nrElems,temp1,temp2,temp3);     
      %TODO here could be an h in RHS
      % [b,temp] = computeRHSwithH(c4n,n4e,s4e,nrSides,area4e, ...
      %   du,tau,Lambda,nrElems,temp1,temp2,temp3,h);     

      %% Solve System
      uNew = zeros(nrSides,1);
      uNew(dof) = A(dof,dof)\b(dof);
      v=(uNew-u)/tau;        

      %% Check Termination
      du = computeGradientNCnew(c4n,n4e,uNew,grads4e);
      ENew = computeEnergy(area4e,uNew,du,alpha,temp,MAMANC);

      dt_u = (u-uNew)/tau; 
      %corr = sqrt(dt_u'*C*dt_u); % Bartels termination criterion
      corr = sqrt(dt_u'*STIMANC*dt_u); % Only gradients
      fprintf('corr/epsStop: %e / %e\n',corr,epsStop);
      format long;
      fprintf('E = %f, E_exact = %f\n', E, -2.05802391003896);
      format short;
      fprintf('============================== \n');

      u = uNew;
      E = ENew;
      energyVec(end+1) = E;
      corrVec(end+1) = corr;

      if saveScreenshots > 0 & mod(length(energyVec),saveScreenshots) == 0
        dirName = sprintf(...
          '../../../../results/tvRegPrimalDualScreenshots/%s/%s',...
          firstScreenshot,datestr(now,'yy_mm_dd_HH_MM_SS'));
        mkdir(dirName);

        name = sprintf('%s/energyVec.txt',dirName);
        file = fopen(name,'w');
        fprintf(file, '%.8g\n',energyVec);
        fclose(file);

        energyRateFig = figure('visible','off');
        loglog(energyVec-energyVec(end));
        ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        title(ftitle);
        xlabel('number of iterations');
        ylabel('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        fName = sprintf('%s/energyRateFigRed%d.png',dirName,red);
        saveas(energyRateFig,fName);

        energyRateFigSemilog = figure('visible','off');
        semilogy(energyVec-energyVec(end));
        ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        title(ftitle);
        xlabel('number of iterations');
        ylabel('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        fName = sprintf('%s/energyRateFigSemilogRed%d.png',dirName,red);
        saveas(energyRateFigSemilog,fName);
      end

      % plotCR(c4n,n4e,uNew);
      % clf('reset');
      % fprintf('\n')
  end
end
