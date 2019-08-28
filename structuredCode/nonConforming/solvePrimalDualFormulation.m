function  [u,corrVec,energyVec] = ...
    solvePrimalDualFormulation(params, currData, u, varLambda) 
% TODO 
% This does ...
%
% solvePrimalDualFormulation.m
% input: 
%
% output:
  
%function  [u,corrVec,energyVec,nrDof] = ...
%  tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,parTau,red,epsStop,parAlpha,f,u,varLambda, ...
%                  saveScreenshots) 

  % extract necessary data
  parTau = params.parTau;
  parAlpha = params.parAlpha;
  showProgress = params.showProgress;

  c4n = currData.c4n;
  n4e = currData.n4e;
  nrSides = currData.nrSides;
  epsStop = currData.epsStop;
  gradsCR4e = currData.gradsCR4e;  
  nrElems = currData.nrElems;
  s4e = currData.s4e;
  area4e = currData.area4e;
  dof = currData.dof;

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

      varLambda = M./repmat(sqrt(sum(M.^2, 2)), 1, 2);
      varLambda(isinf(varLambda)) = 0;
      
      % compute RHS
      
      % TODO recalculate this at some point, Tiens looks different (smarter prob.)

      b = zeros(nrSides, 1);
      rhsInt = zeros(nrSides, 1); % int_\Omega f*u dx
      for elem = 1 : nrElems
        bLocal = (gradCRu(elem, :)/parTau - ...
          varLambda(elem, :))*gradsCR4e(:, :, elem)';  
        rhsInt(s4e(elem, :)) = rhsInt(s4e(elem,:)) + ...
          [rhsInt(elem), rhsInt(elem), rhsInt(elem)]';
        b(s4e(elem, :)) = b(s4e(elem, :)) + area4e(elem)*bLocal'; % right-hand side
      end
      b = b + rhsInt;

      %% Solve System
      uNew = zeros(nrSides, 1);
      uNew(dof) = A(dof, dof)\b(dof);
      v = (uNew - u)/parTau;        

      %% Check Termination
      gradCRu = gradientCR(currData, u);
      % TODO continue here
      ENew = computeEnergy(area4e,uNew,gradCRu,parAlpha,rhsInt,maMaCR);

      % TODO here case distinction for termination criteria needs to be done
      dt_u = (u - uNew)/parTau; 
      %corr = sqrt(dt_u'*C*dt_u); % Bartels termination criterion
      % TODO find out what the hell corr is, how does Tien call this
      corr = sqrt(dt_u'*stiMaCR*dt_u); % Only gradients

      if showProgress
        fprintf('corr/epsStop: %e / %e\n',corr,epsStop);
        format long;
        fprintf('E = %f, E_exact = %f\n', E, -2.05802391003896);
        format short;
        fprintf('============================== \n');
      end

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
