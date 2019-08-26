%TODO what is the name of this
function  [u,corrVec,energyVec] = ...
    tvRegPrimalDual(params, currData, u, varLambda) 
% TODO 
%
% tvRegPrimalDual.m
% input: 
%
% output:
%        
  
%function  [u,corrVec,energyVec,nrDof] = ...
%  tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,parTau,red,epsStop,parAlpha,f,u,varLambda, ...
%                  saveScreenshots) 

  % extract necessary data
  parTau = params.parTau;
  parAlpha = params.parAlpha;
  c4n = currData.c4n;
  n4e = currData.n4e;
  nrSides = currData.nrSides;
  epsStop = currData.epsStop;
  gradsCR4e = currData.gradsCR4e;  

  % INITIALIZATION 
  
  % extract necessary data

  % initialize remaing parameters
  firstScreenshot = datestr(now, 'yy_mm_dd_HH_MM_SS');

  [stiMaCR, maMaCR] = computeFeMatricesCR(currData);

  A = stiMaCR/parTau+parAlpha*maMaCR; %TODO here could be an h in front of stiMaCR
  %C = maMaCR + h*stiMaCR;

  gradCRu = gradientCR(currData, u);

  v = zeros(nrSides, 1);    

  corr = epsStop+1; 
  corrVec = [];
  energyVec = [];
  E = 1;

  while corr > epsStop
      gradCRv = gradientCR(currData, v);
      M = varLambda + parTau*(gradCRu + parTau*gradCRv);

      varLambda = M./repmat(sqrt(sum(M.^2,2)),1,2);
      varLambda(isinf(varLambda)) = 0;
      
      % TODO continue here
      [b,rhsInt] = computeRHS(c4n,n4e,s4e,nrSides,area4e, ...
        gradCRu,parTau,varLambda,nrElems,rhsInt1,rhsInt2,rhsInt3);     
      %TODO here could be an h in RHS
      % [b,temp] = computeRHSwithH(c4n,n4e,s4e,nrSides,area4e, ...
      %   du,parTau,varLambda,nrElems,temp1,temp2,temp3,h);     

      %% Solve System
      uNew = zeros(nrSides,1);
      uNew(dof) = A(dof,dof)\b(dof);
      v=(uNew-u)/parTau;        

      %% Check Termination
      gradCRu = gradientCR(currData, u);
      ENew = computeEnergy(area4e,uNew,gradCRu,parAlpha,rhsInt,maMaNC);

      dt_u = (u-uNew)/parTau; 
      %corr = sqrt(dt_u'*C*dt_u); % Bartels termination criterion
      corr = sqrt(dt_u'*stiMaCR*dt_u); % Only gradients
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
%TODO write function to save all the stuff from the screenshot, just so this 
%code gets only half as long (and gets more readable + gt)
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

      if showPlots
        plotCR(c4n,n4e,uNew);
        clf('reset');
        fprintf('\n')
      end
  end
end
