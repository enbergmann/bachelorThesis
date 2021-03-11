function plotElasticityLS(configName, problemData, geometryData, levelData, parameter, figStruct, solutionData, exSolutionData, STIMA, fMean4e)
  
  plotLable = configName;
  plotLable = strrep(plotLable,'c_','');
  plotLable = strrep(plotLable,'default','');
  plotLable = strrep(plotLable,'_',' und ');
  plotLable = strrep(plotLable,'^','\^');
  
  c4n = geometryData.c4n;
  n4e = geometryData.n4e;
  nrDof4lvl = levelData.nrDof4lvl;
  
  n4s      = computeN4s(n4e);
  nrNodes  = size(c4n,1);
  nrSides  = size(n4s,1);
  mid4e    = computeMid4e(c4n, n4e);
  s4e      = computeS4e(n4e);
  mid4e    = computeMid4e(c4n,n4e);
  signs4e  = computeSigns4e(c4n,n4e,n4s,s4e);
  length4s = computeLength4s(c4n,n4s);
  coord4e  = [c4n(n4e(:,1),:),c4n(n4e(:,2),:),c4n(n4e(:,3),:)];
  area4e   = computeArea4e(c4n,n4e);
  mid4s    = computeMid4s(c4n, n4s);
  length4e = length4s(s4e);
  %% convergence plot
  if parameter.convergencePlot
    figStruct.convFig = FigureAdvanced(figStruct.convFig);
    if ~strcmp(parameter.solver,'P1P1')
      plotConvergence(nrDof4lvl,levelData.residual4lvl,'Residual');
      hold all;
    end
    if parameter.showLSdetail
      if ~strcmp(parameter.solver,'P1P1')
        plotConvergence(nrDof4lvl,levelData.contrDiv4lvl,'||f+div(\sigma_{h})||');
      end
      plotConvergence(nrDof4lvl,levelData.contrC4lvl,'||C^{s}\sigma_{h}-C^{1+s}\epsilon(u_{h})||');
      hold all
    end
    if parameter.showQ
      plotConvergence(nrDof4lvl,levelData.q4lvl,'q');
      hold all
    end
    if parameter.computeDataError && parameter.showDataError && levelData.dataApprox4lvl(1) > 10^-8
        plotConvergence(nrDof4lvl,levelData.dataApprox4lvl,'||f+f_l||');
        hold all
    end
    if parameter.computeQuasiExactError
      plotConvergence(nrDof4lvl,levelData.quasiExactErr4lvl,'quasi exact error');
      hold all
    end
    if problemData.exSolKnown && parameter.useExactSolution
      if parameter.showErrorDetail
        plotConvergence(nrDof4lvl,levelData.errorU4lvl,'||u-u_{h}||_{H^1}');
        hold all
        plotConvergence(nrDof4lvl,levelData.errorL2U4lvl,'||u-u_{h}||');
        plotConvergence(nrDof4lvl,levelData.errorDu4lvl,'||Du-Du_{h}||');
        plotConvergence(nrDof4lvl,levelData.errorL2Sigma4lvl,'||\sigma-\sigma_{h}||');
        if strcmp(parameter.solver,'P1P1')
          plotConvergence(nrDof4lvl,levelData.errorSigmaP04lvl,'||\sigma-\sigma_0||');
        end
      end
      if ~strcmp(parameter.solver,'P1P1')
        plotConvergence(nrDof4lvl,levelData.error4lvl,'exact error');
        hold all
        if parameter.showErrorDetail
          plotConvergence(nrDof4lvl,levelData.errorSigma4lvl,'||\sigma-\sigma_{h}||_{H(div)}');
          plotConvergence(nrDof4lvl,levelData.errorDivSigma4lvl,'||div(\sigma-\sigma_{h})||');
        end
      end
      if parameter.showElastError
        plotConvergence(nrDof4lvl,levelData.elastError4lvl,'|||u-u_{h}|||');
      end
      if parameter.computeDivU && parameter.showDivU
        plotConvergence(nrDof4lvl,levelData.errorDivU4lvl,'||div (u-u_h)||');
      end
    end
    title({plotLable,'','convergence plot'})
    hold off
  end
  %% mesh plot
  if parameter.meshPlot
    figStruct.meshFig = FigureAdvanced(figStruct.meshFig);
    plotTriangulation(c4n,geometryData.n4e);
  end
  %% solution plot u
  if parameter.solutionPlotU
    figStruct.solutionFigU = FigureAdvanced(figStruct.solutionFigU);
    %% draw exakt solution of u
    if problemData.exSolKnown && parameter.useExactSolution
      nrSubplotsPerRow = 2;
      subplot(2,nrSubplotsPerRow,3); plotP1(c4n,n4e,exSolutionData.u1Cex(c4n)); colorbar; title('exact u first component');
      subplot(2,nrSubplotsPerRow,4); plotP1(c4n,n4e,exSolutionData.u2Cex(c4n)); colorbar; title('exact u second component');
    else
      nrSubplotsPerRow = 1;
    end
    %% draw FEM solution of u
    switch parameter.method
      case 'S1'
        subplot(2,nrSubplotsPerRow,1); plotP1(c4n,n4e,solutionData.u(1:2:end)); colorbar; title({plotLable,'','u first component'});
        subplot(2,nrSubplotsPerRow,2); plotP1(c4n,n4e,solutionData.u(2:2:end)); colorbar; title('u second component');
      case 'KS'
        subplot(2,nrSubplotsPerRow,1); plotP1(c4n,n4e,solutionData.u(1:nrNodes)); colorbar; title({plotLable,'','u first component'});
        subplot(2,nrSubplotsPerRow,2); plotCR(c4n,n4e,solutionData.u(nrNodes+1:nrNodes + nrSides)); colorbar; title('u second component');
      case 'CR'
        subplot(2,nrSubplotsPerRow,1); plotCR(c4n,n4e,solutionData.u(1:2:end)); colorbar; title('u first component');
        subplot(2,nrSubplotsPerRow,2); plotCR(c4n,n4e,solutionData.u(2:2:end)); colorbar; title('u second component');
    end
   end
  %% solution plot sigma
  if parameter.solutionPlotSigma
    figStruct.solutionFigSigma = FigureAdvanced(figStruct.solutionFigSigma);
    %% draw exakt solution of sigma
    if problemData.exSolKnown && parameter.useExactSolution
      sigmaEx4mid = exSolutionData.sigmaExact(mid4e);
      nrSubplotsPerRow = 2;
      subplot(2,nrSubplotsPerRow,3);
        quiver2(mid4e(:,1),mid4e(:,2),reshape(sigmaEx4mid(1,1,:),size(n4e,1),1),reshape(sigmaEx4mid(1,2,:),size(n4e,1),1),'n=',0.1,'w=',[1 1]);
        colorbar; title('exact sigma first component');
      subplot(2,nrSubplotsPerRow,4);
        quiver2(mid4e(:,1),mid4e(:,2),reshape(sigmaEx4mid(2,1,:),size(n4e,1),1),reshape(sigmaEx4mid(2,2,:),size(n4e,1),1),'n=',0.1,'w=',[1 1]);
        colorbar; title('exact sigma second component');
    else
      nrSubplotsPerRow = 1;
    end
    %% draw FEM solution of sigma
    sigma4e   = solutionData.sigma4e; %solutionData{10};
    sigmaRT4e = ...
                [-sum(signs4e.*length4e.*sigma4e(:,1:2:5).*coord4e(:,[5 1 3]),2),...
                 -sum(signs4e.*length4e.*sigma4e(:,1:2:5).*coord4e(:,[6 2 4]),2),...
                  sum(signs4e.*length4e.*sigma4e(:,1:2:5),2),...
                 -sum(signs4e.*length4e.*sigma4e(:,2:2:6).*coord4e(:,[5 1 3]),2),...
                 -sum(signs4e.*length4e.*sigma4e(:,2:2:6).*coord4e(:,[6 2 4]),2),...
                  sum(signs4e.*length4e.*sigma4e(:,2:2:6),2)] ./ repmat(2*area4e,1,6);
    subplot(2,nrSubplotsPerRow,1);
      quiver2(mid4e(:,1),mid4e(:,2),sigmaRT4e(:,1)+sigmaRT4e(:,3).*mid4e(:,1), sigmaRT4e(:,2)+sigmaRT4e(:,3).*mid4e(:,2),'n=',0.1,'w=',[1 1]);
      colorbar; title({plotLable,'','sigma first component'});
    subplot(2,nrSubplotsPerRow,2);
      quiver2(mid4e(:,1),mid4e(:,2),sigmaRT4e(:,4)+sigmaRT4e(:,6).*mid4e(:,1), sigmaRT4e(:,5)+sigmaRT4e(:,6).*mid4e(:,2),'n=',0.1,'w=',[1 1]);
      colorbar; title('sigma second component');
  end
  %% displacement plot
  if parameter.displacementPlot
    figStruct.displacementFig = FigureAdvanced(figStruct.displacementFig);
    %% draw exakt displacement
    if problemData.exSolKnown && parameter.useExactSolution
      nrSubplotsPerColumn = 2;
      subplot(nrSubplotsPerColumn,1,2);
        quiver2(mid4e(:,1),mid4e(:,2),exSolutionData.u1Cex([mid4e(:,1),mid4e(:,2)]),exSolutionData.u2Cex([mid4e(:,1),mid4e(:,2)]),'n=',0.1,'w=',[1 1]);
        colorbar; title({plotLable,'','exact displacement'});
    else
      nrSubplotsPerColumn = 1;
    end
    %% draw FEM displacement
    switch parameter.method
      case 'S1'
        subplot(nrSubplotsPerColumn,1,1);
          quiver2(c4n(:,1),c4n(:,2),solutionData.u(1:2:end),solutionData.u(2:2:end),'n=',0.1,'w=',[1 1]);
          title({plotLable,'','displacement'});
      case 'KS'
        subplot(nrSubplotsPerColumn,1,1);
          u1CatMid4s = (solutionData.u(n4s(:,1))+solutionData.u(n4s(:,2)))/2;
          quiver2(mid4s(:,1),mid4s(:,2),u1CatMid4s(:),solutionData.u(1+nrNodes:nrNodes+nrSides),'n=',0.1,'w=',[1 1]);
          title({plotLable,'','displacement'});
      case 'CR'
        subplot(nrSubplotsPerColumn,1,1);
          title({plotLable,'','displacement (missing)'});
    end
  end
  %% deformation plot
  if parameter.deformationPlot
    figStruct.deformationFig = FigureAdvanced(figStruct.deformationFig);
    switch parameter.method
      case 'S1'
        plotP1P1(n4e,c4n,solutionData.u,parameter.lambda,parameter.lame_mu,1);
      case 'KS'
        title('missing');
    end
  end
  %% STIMA plot
  if parameter.STIMAPlot
    figStruct.STIMAFig = FigureAdvanced(figStruct.STIMAFig);
    spy(STIMA);
  end
  %% time-residual plot
  if parameter.timePlot
    figStruct.timeFig = FigureAdvanced(figStruct.timeFig);
    plotConvergence(levelData.time4lvl,levelData.residual4lvl,'Residual');
    hold all;
    plotConvergence(levelData.time4lvl,levelData.contrDiv4lvl,'||f+div(\sigma_{h})||');
    plotConvergence(levelData.time4lvl,levelData.contrC4lvl,'||C^{s}\sigma_{h}-\epsilon(u_{h})||');
  end
  %% rhs plot
  if parameter.rhsPlot && parameter.computeDataError
    figStruct.rhsFig = FigureAdvanced(figStruct.rhsFig);
    N = 200;
    [X, Y]   = meshgrid(linspace(min(min(c4n)),max(max(c4n)), N));
    Z = problemData.f([X(:), Y(:)]);
    subplot(2,2,1); surf(X, Y, reshape(Z(:, 1), N, N),'LineStyle','none'); title('exact f first component');
    subplot(2,2,2); surf(X, Y, reshape(Z(:, 2), N, N),'LineStyle','none'); title('exact f second component');
    subplot(2,2,3); plotP04e(c4n,n4e,fMean4e(:,1)); title('P_0 approx of f first component');
    subplot(2,2,4); plotP04e(c4n,n4e,fMean4e(:,2)); title('P_0 approx of f second component');
  end
  %% condition plot
  if parameter.condPlot && parameter.computeCond
    figStruct.condFig = FigureAdvanced(figStruct.condFig);
    plotConvergence(nrDof4lvl,levelData.cond4lvl,'condition');
  end
  %% divergence u plot
  if parameter.divUPlot && parameter.computeDivU
  % draw exact divergence of u
    figStruct.divFig = FigureAdvanced(figStruct.divFig);
    if problemData.exSolKnown && parameter.useExactSolution
      nrSubplotsPerColumn = 2;
      divUexact = exSolutionData.divUexact; %exSolutionData{10};
      Du1CexDx  = exSolutionData.Du1CexDx; %exSolutionData{3};
      Du2CexDy  = exSolutionData.Du2CexDy; %exSolutionData{6};
      subplot(nrSubplotsPerColumn,3,4); plotP04e(c4n,n4e,exSolutionData.divUexact(mid4e)); title('exact divergence of u');
      subplot(nrSubplotsPerColumn,3,5); plotP04e(c4n,n4e,exSolutionData.Du1CexDx(mid4e)); title('exact partial derivative d_{x_1}u_1');
      subplot(nrSubplotsPerColumn,3,6); plotP04e(c4n,n4e,exSolutionData.Du2CexDy(mid4e)); title('exact partial derivative d_{x_2}u_2');
    else
      nrSubplotsPerColumn = 1;
    end
  % draw FEM divergence of u
    divU4e   = solutionData.divU4e; %solutionData{8};
    Du1CDx4e = solutionData.Du1CDx4e; %solutionData{2};
    Du2CDy4e = solutionData.Du2CDy4e; %solutionData{5};
    subplot(nrSubplotsPerColumn,3,1); plotP04e(c4n,n4e,solutionData.divU4e); title({plotLable,'','discrete divergence of u'});
    subplot(nrSubplotsPerColumn,3,2); plotP04e(c4n,n4e,solutionData.Du1CDx4e); title('discrete partial derivative d_{x_1}u_1');
    subplot(nrSubplotsPerColumn,3,3); plotP04e(c4n,n4e,solutionData.Du2CDy4e); title('discrete partial derivative d_{x_2}u_2');
  end
  %% divergence sigma plot
  if parameter.divSigmaPlot
    figStruct.divSigmaFig = FigureAdvanced(figStruct.divSigmaFig);
    if problemData.exSolKnown && parameter.useExactSolution
      nrSubplotsPerColumn = 2;
      fval = problemData.f(mid4e);
      f1Cval = fval(:,1); f2Cval = fval(:,2);
      subplot(nrSubplotsPerColumn,3,4); plotP04e(c4n,n4e,-f1Cval);
        title('exact divergence of sigma (first component) (-f_1)')
      subplot(nrSubplotsPerColumn,3,5); plotP04e(c4n,n4e,-f2Cval);
        title('exact divergence of sigma (second component) (-f_2)')
      subplot(nrSubplotsPerColumn,3,6);  
        quiver2(mid4e(:,1),mid4e(:,2),-f1Cval,-f2Cval,'n=',0.1,'w=',[1 1]);
        title('exact divergence of sigma (-f)');
    else
      nrSubplotsPerColumn = 1;
    end
    subplot(nrSubplotsPerColumn,3,1); plotP04e(c4n,n4e,solutionData.divSigma4e(:,1));
      title({plotLable,'','discrete divergence of \sigma (first component)'});
    subplot(nrSubplotsPerColumn,3,2); plotP04e(c4n,n4e,solutionData.divSigma4e(:,2));
      title('discrete divergence of sigma (second component)');
    subplot(nrSubplotsPerColumn,3,3);
      quiver2(mid4e(:,1),mid4e(:,2),solutionData.divSigma4e(:,1),solutionData.divSigma4e(:,2),'n=',0.1,'w=',[1 1]);
      title('discrete divergence of sigma');
  end
  %% overlay plot
  if parameter.overlayPlot && parameter.computeQuasiExactError
     %figStruct.overlayFig = FigureAdvanced(figStruct.overlayFig);
    %subplot(3,4,1);  plotTriangulation(c4n,n4e); title({plotLable,'','mesh 1'});
    %subplot(3,4,3);  plotTriangulation(c4nF,n4eF); title('mesh 2');
    %subplot(3,4,4);  plotTriangulation(c4n_overlay,n4e_overlay); title('mesh overlay');
    %subplot(3,4,5);  plotP1(c4n,n4e,u(1:2:end)); title('u1C');
    %subplot(3,4,6);  plotP1(c4n,n4e,u(2:2:end)); title('u2C');
    %subplot(3,4,7);  plotP1(c4nF,n4eF,uF(1:2:end)); title('u_{qex}1C');
    %subplot(3,4,8);  plotP1(c4nF,n4eF,uF(2:2:end)); title('u_{qex}2C');
    %subplot(3,4,9);  plotP1(c4n_overlay,n4e_overlay,u1C_overlay); title('u1C\_overlay');
    %subplot(3,4,10); plotP1(c4n_overlay,n4e_overlay,u2C_overlay); title('u2C\_overlay');
    %subplot(3,4,11); plotP1(c4n_overlay,n4e_overlay,uF1C_overlay); title('u_{qex}1C\_overlay');
    %subplot(3,4,12); plotP1(c4n_overlay,n4e_overlay,uF2C_overlay); title('u_{qex}1C\_overlay');
  end
  drawnow
end
