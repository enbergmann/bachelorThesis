function  [u,corrVec,energyVec,nrDoF] = ...
    tvRegPrimalDual(params, c4n, n4e, n4sDb, n4sNb, u, varLambda,...
    epsStop, h, red) 

  
%function  [u,corrVec,energyVec,nrDoF] = ...
%  tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,parTau,red,epsStop,parAlpha,f,u,varLambda, ...
%                  saveScreenshots) 

  %saveScreenshots: integer x bigger than zero (save results every x iterations)
  %                 or equal to zero (don't save screenshots)

  % unpack params
  parTau = params.parTau;
  parAlpha = params.parAlpha;
  f = params.f;
  saveScreenshots = params.saveScreenshots;
 
  
  firstScreenshot = datestr(now,'yy_mm_dd_HH_MM_SS');

  nrElems = size(n4e,1);
  area4e = computeArea4e(c4n,n4e);
  s4e = computeS4e(n4e);
  nrSides = max(max(s4e));

  dof = computeDof(n4e,nrSides,n4sDb,n4sNb);
  nrDoF = length(dof);

  [stiMaNC,maMaNC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);
  A = stiMaNC/parTau+parAlpha*maMaNC; %TODO here could be an h in front of stiMaNC
  C = maMaNC + h*stiMaNC;

  [rhsInt1,rhsInt2,rhsInt3] = computeIntegrals(f,c4n,n4e,200,area4e);

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
      gradsT = [ones(1,3);c4n(n4e(elem,:),:)']\[zeros(1,2);-2*eye(2)];
      gradsT = gradsT([3 1 2],:);
      grads4e(:,:,elem) = gradsT;
  end



  
  while corr > epsStop
      dv = computeGradientNCnew(c4n,n4e,v,grads4e);
      M = varLambda + parTau*(du + parTau*dv);
      varLambda = bsxfun(@rdivide,M,max(1,sqrt(sum(M.^2,2))));
      
      [b,rhsInt] = computeRHS(c4n,n4e,s4e,nrSides,area4e, ...
        du,parTau,varLambda,nrElems,rhsInt1,rhsInt2,rhsInt3);     
      %TODO here could be an h in RHS
      % [b,temp] = computeRHSwithH(c4n,n4e,s4e,nrSides,area4e, ...
      %   du,parTau,varLambda,nrElems,temp1,temp2,temp3,h);     

      %% Solve System
      uNew = zeros(nrSides,1);
      uNew(dof) = A(dof,dof)\b(dof);
      v=(uNew-u)/parTau;        

      %% Check Termination
      du = computeGradientNCnew(c4n,n4e,uNew,grads4e);
      ENew = computeEnergy(area4e,uNew,du,parAlpha,rhsInt,maMaNC);

      dt_u = (u-uNew)/parTau; 
      %corr = sqrt(dt_u'*C*dt_u); % Bartels termination criterion
      corr = sqrt(dt_u'*stiMaNC*dt_u); % Only gradients
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

        limitAitken = ...
          (energyVec(end)*energyVec(end-2)-energyVec(end-1)^2)/...
          (energyVec(end-2)-2*energyVec(end-1)+energyVec(end));

        name = sprintf('%s/AitkenLimitEnergyVec.txt',dirName);
        file = fopen(name,'w');
        fprintf(file, 'Aitken Limit\n%.8g\n',limitAitken);
        fprintf(file, '%.8g\n',energyVec);
        fclose(file);

        % Use Aitken Extrapolation as limit
        energyRateAitkenFig = figure('visible','off');
        loglog(energyVec-limitAitken);
        ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
        title(ftitle);
        xlabel('number of iterations');
        ylabel('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
        fName = sprintf('%s/EnergyRateFigRed%dAitken.png',dirName,red);
        saveas(energyRateAitkenFig,fName);

        energyRateAitkenFigSemilog = figure('visible','off');
        semilogy(energyVec-limitAitken);
        ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
        title(ftitle);
        xlabel('number of iterations');
        ylabel('|E_{NC}(u_{j})-E_{NC}(u_{aitken})|');
        fName = sprintf('%s/EnergyRateFigSemilogRed%dAitken.png',dirName,red);
        saveas(energyRateAitkenFigSemilog,fName);

        % Use last energyVec value as limit
        energyRateFig = figure('visible','off');
        loglog(energyVec-energyVec(end));
        ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        title(ftitle);
        xlabel('number of iterations');
        ylabel('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        fName = sprintf('%s/EnergyRateFigRed%d.png',dirName,red);
        saveas(energyRateFig,fName);

        energyRateFigSemilog = figure('visible','off');
        semilogy(energyVec-energyVec(end));
        ftitle=sprintf('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        title(ftitle);
        xlabel('number of iterations');
        ylabel('|E_{NC}(u_{j})-E_{NC}(u_{end})|');
        fName = sprintf('%s/EnergyRateFigSemilogRed%d.png',dirName,red);
        saveas(energyRateFigSemilog,fName);
      end

      % plotCR(c4n,n4e,uNew);
      % clf('reset');
      % fprintf('\n')
  end
end
