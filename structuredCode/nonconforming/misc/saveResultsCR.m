function saveResultsCR(params, currData, ...
    outputLvlInfo, outputLvlError, outputLvlEnergy, outputLvlHidden, output)

%% INIT
  % extract necessary parameters from params
  geometry = params.geometry;
  f = params.f;
  benchmark = params.benchmark;
  exactSolutionKnown = params.exactSolutionKnown;
  uExact = params.uExact;
  expName = params.expName;
  dirInfoName = params.dirInfoName;
  figVisible = params.figVisible;
  parAlpha = params.parAlpha;
  parBeta = params.parBeta;
  rhsGradientKnown = params.rhsGradientKnown;
  useExactEnergy = params.useExactEnergy;
  exactEnergy = params.exactEnergy;
  plotGivenFunctions = params.plotGivenFunctions;
  refinementLevel4Plots = params.refinementLevel4Plots;
  polygonMesh = params.polygonMesh;
  parTheta = params.parTheta;
  minNrDof = params.minNrDof;
  d = params.d;
  parGamma = params.parGamma;

  % extract necessary information from currData
  nrElems = currData.nrElems; 
  nrSides = currData.nrSides; 
  nrDof = currData.nrDof; 
  hMax = currData.hMax;
  epsStop = currData.epsStop;
  c4n = currData.c4n;
  n4e = currData.n4e;
  s4e = currData.s4e;
  n4sDb = currData.n4sDb;
  n4sNb = currData.n4sNb;

  % extract necessary information from outputLvlInfo
  currLvl = outputLvlInfo.lvl(end);
  nrDof4lvl = outputLvlInfo.nrDof;
  time4lvl = outputLvlInfo.time;
  time = time4lvl(end);
  nrIterations4lvl = outputLvlInfo.nrIterations;

  % extract necessary information from outputLvlError
  eta4lvl = outputLvlError.eta;
  etaVol4lvl = outputLvlError.etaVol;
  etaJumps4lvl = outputLvlError.etaJumps;
  if exactSolutionKnown, error4lvl = outputLvlError.error4lvl; end

  % extract necessary information from outputLvlInfoEnergy
  if useExactEnergy, diffDiscExactE4lvl = outputLvlEnergy.diffDiscExactE; end
  if rhsGradientKnown 
    diffGuebGleb4lvl = outputLvlEnergy.diffGuebGleb; 
    if useExactEnergy, diffExactEGleb4lvl = outputLvlEnergy.diffExactEGleb; end
  end
  
  % extract necessary information from output
  u = output.u;
  corrVec = output.corrVec;
  energyVec = output.energyVec;
  eNcAbsDiffVec = output.otherCorr.eNcAbsDiffVec;
  bar15TerminationVec = output.otherCorr.bar15TerminationVec;
  bar15TerminationWithoutL2Vec = output.otherCorr.bar15TerminationWithoutL2Vec; 
  bar12TerminationSqrtVec = output.otherCorr.bar12TerminationSqrtVec; 
   
%% MAIN
  % create directory
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf('../../results/nonconforming/%s/%s/lvl%d', ...
    expName, dirInfoName, currLvl);
    % startAlgorithmNC run from .../structuredcode/nonconforming
  if nrDof<1e6 % to avoid clogging results with large files
    mkdir(sprintf('%s/mesh', dirName)); 
  end
  mkdir(sprintf('%s/iteration', dirName));
  warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
  % save information about the experiment
  if currLvl == 0
    % this means the benchmark-file should not be changed until level 0 is
    % saved
    source = sprintf('benchmarks/%s.m', benchmark);
    destination = ...
      sprintf('../../results/nonconforming/%s/%s/benchmark_%s.m', ...
      expName, dirInfoName, benchmark);
    copyfile(source, destination);

    % save experiment parameters
    expParams = struct(...
      'parAlpha', parAlpha, 'parBeta', parBeta, ...
      'minNrDof', minNrDof, 'parTheta', parTheta, ...
      'd', d, 'parGamma', parGamma);
    name = sprintf('../../results/nonconforming/%s/%s/expParams.csv', ...
      expName, dirInfoName);
    writetable(struct2table(expParams, 'AsArray', true), name);

    % save remaining parameters
    paramsReduced = rmfield(params, ...
      {'c4n', 'n4e', 'n4sDb', 'n4sNb', 'parTheta', 'minNrDof', ...
      'd', 'parGamma', 'parAlpha', 'parBeta'});
    name = sprintf('../../results/nonconforming/%s/%s/paramsReduced.csv', ...
      expName, dirInfoName);
    writetable(struct2table(paramsReduced, 'AsArray', true), name);

    % plot rhs and grayscale image of rhs
    if plotGivenFunctions
      if polygonMesh
        [c4nRhs, n4eRhs] = computeGeometryPolygon(refinementLevel4Plots);
      else
        [c4nRhs, n4eRhs] = loadGeometry(geometry, refinementLevel4Plots);
      end

      fVal = f(c4nRhs);

      rhsFig = figure('visible', figVisible); 
      trisurf(n4eRhs, c4nRhs(:, 1), c4nRhs(:, 2), fVal, ...
        'EdgeColor', 'None');
      fName = sprintf('../../results/nonconforming/%s/%s/rhs.png', ...
        expName, dirInfoName);
      saveas(rhsFig, fName);

      rhsAxisFig = figure('visible', figVisible); 
      plotAxis(c4nRhs, fVal);
      fName = sprintf('../../results/nonconforming/%s/%s/rhsAxis.png', ...
        expName, dirInfoName);
      saveas(rhsAxisFig, fName);

      rhsGrayscaleFig = figure('visible', figVisible); 
      trisurf(n4eRhs, c4nRhs(:, 1), c4nRhs(:, 2),  fVal, 'EdgeColor', 'None');
      view(0, 90);
      axis image;
      colormap gray;
      fName = sprintf('../../results/nonconforming/%s/%s/rhsGrayscale.png', ...
        expName, dirInfoName);
      saveas(rhsGrayscaleFig, fName);

      if exactSolutionKnown
        uExactVal = uExact(c4nRhs);

        uExactFig = figure('visible', figVisible); 
        trisurf(n4eRhs, c4nRhs(:, 1), c4nRhs(:, 2), uExactVal, ...
          'EdgeColor', 'None');
          sprintf('../../results/nonconforming/%s/%s/exactSolution.png', ...
          expName, dirInfoName);
        fName = sprintf(...
          '../../results/nonconforming/%s/%s/exactSolution.png', ...
          expName, dirInfoName);
        saveas(uExactFig, fName);


        uExactAxisFig = figure('visible', figVisible); 
        plotAxis(c4nRhs, uExactVal);
        fName = sprintf(...
          '../../results/nonconforming/%s/%s/exactSolutionAxis.png', ...
          expName, dirInfoName);
        saveas(uExactAxisFig, fName);

        uExactGrayscaleFig = figure('visible', figVisible); 
        trisurf(n4eRhs, c4nRhs(:, 1), c4nRhs(:, 2), uExactVal, ...
          'EdgeColor', 'None');
        view(0, 90);
        axis image;
        colormap gray;
        fName = sprintf(...
          '../../results/nonconforming/%s/%s/exactSolutionGrayscale.png', ...
          expName, dirInfoName);
        saveas(uExactGrayscaleFig, fName);
      end
    end
  end

  % save parameters for current level
  currDataReduced = struct(...
    'nrElems', nrElems, 'nrNodes', size(c4n, 1), 'nrSides', nrSides, ...
    'nrDof', nrDof, 'hMax', hMax, 'epsStop', epsStop);
  name = sprintf('%s/currentDataReduced.csv', dirName);
  writetable(struct2table(currDataReduced, 'AsArray', true), name);

  % save plots of solution
  approxFig = figure('visible', figVisible); 
  plotCR(c4n, n4e, u);
  fName = sprintf('%s/solution.png', dirName);
  saveas(approxFig, fName);
  
  approxFigAxis = figure('visible', figVisible); 
  plotAxisNC(c4n,n4e,u);
  fName = sprintf('%s/solutionAxis.png', dirName);
  saveas(approxFigAxis, fName);

  grayscaleFig = figure('visible', figVisible); 
  plotGrayscale(c4n, n4e, mean(u(s4e), 2));
  fName = sprintf('%s/solutionGrayscale.png', dirName);
  saveas(grayscaleFig, fName);

  % save plots and results of the iteration for the level
  % save corrections
  name = sprintf('%s/iteration/corrVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', corrVec);
  fclose(file);
  
  corrFig = figure('visible', figVisible);
  loglog(corrVec);
  ftitle = sprintf(...
    'loglog plot - corr for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('corr');
  fName = sprintf('%s/iteration/corr.png', dirName);
  saveas(corrFig, fName);

  name = sprintf('%s/iteration/eNcAbsDiffVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', eNcAbsDiffVec);
  fclose(file);

  name = sprintf('%s/iteration/bar15TerminationVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', bar15TerminationVec);
  fclose(file);

  name = sprintf('%s/iteration/bar15TerminationWithoutL2Vec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', bar15TerminationWithoutL2Vec);
  fclose(file);

  name = sprintf('%s/iteration/bar12TerminationSqrtVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', bar12TerminationSqrtVec);
  fclose(file);

  terminationFig = figure('visible', figVisible);
  loglog(corrVec);
  hold on
  loglog(eNcAbsDiffVec);
  loglog(bar15TerminationVec);
  loglog(bar15TerminationWithoutL2Vec);
  loglog(bar12TerminationSqrtVec);
  hold off
  ftitle = sprintf(...
    'loglog plot - corr for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  xlabel('number of iterations');
  legend([sprintf("corr (energyNormDiffernce)"), ...
    sprintf("eNcAbsDiff"), sprintf("bar15"), ...
    sprintf("bar15TerminationWithoutL2"), sprintf("bar12sqrt")], ...
    'Location', 'SW');
  fName = sprintf('%s/iteration/termination.png', dirName);
  saveas(terminationFig, fName);
  
  % save discrete energies
  name = sprintf('%s/iteration/energyVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', energyVec);
  fclose(file);
  
  enFig = figure('visible', figVisible); 
  enFigLegend = sprintf("nrDof = %d (%0.2fs)", nrDof, time);
  plot(energyVec);
  if useExactEnergy
    hold on;
    plot(exactEnergy*ones(1, length(energyVec)));
    enFigLegend(end+1) = sprintf("E_u = %.8g", exactEnergy);
  end
  legend(enFigLegend);
  ftitle=sprintf('Energy for nrDof=%d, \\alpha =%d, \\beta =%d',...
  nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/iteration/energy.png', dirName);
  xlabel('number of iterations');
  ylabel('energy');
  saveas(enFig, fName);

  if useExactEnergy
    % save differences between discrete energies and exact energy
    name = sprintf('%s/iteration/exactAbsoluteEnergyDifference.txt', dirName);
    file = fopen(name, 'w');
    fprintf(file, '%.8g\n', abs(energyVec-exactEnergy));
    fclose(file);

    enDiffExactFig = figure('visible',figVisible);
    loglog(abs(energyVec-exactEnergy));
    ftitle = sprintf(...
      '|E_{NC}(u_{NC})-E_u| for nrDof = %d, \\alpha = %d, \\beta = %d', ...
      nrDof, parAlpha, parBeta);
    title(ftitle);
    xlabel('number of iterations');
    ylabel('|E_{NC}(u_{NC})-E_u|');
    fName = sprintf('%s/iteration/enDiffExact.png', dirName);
    saveas(enDiffExactFig, fName);
  end

  % save AFEM results and triangulation
  % plot triangulation
  triangFig = figure('visible',figVisible);
  plotTriangulation(c4n,n4e);
  fName = sprintf('%s/triangulation.png', dirName);
  saveas(triangFig, fName);

  if nrDof<1e6 % to avoid clogging results with large files
    dlmwrite(sprintf('%s/mesh/c4n.txt', dirName), c4n, 'Delimiter', '\t');
    dlmwrite(sprintf('%s/mesh/n4e.txt', dirName), n4e, 'Delimiter', '\t');
    dlmwrite(sprintf('%s/mesh/n4eM.txt', dirName), n4e-1, 'Delimiter', '\t');
    dlmwrite(sprintf('%s/mesh/n4sDb.txt', dirName), n4sDb, 'Delimiter', '\t');
    dlmwrite(sprintf('%s/mesh/n4sNb.txt', dirName), n4sNb, 'Delimiter', '\t');
    dlmwrite(sprintf('%s/mesh/s4e.txt', dirName), s4e, 'Delimiter', '\t');
  end

  % save AFEM results
  tableStruct = outputLvlInfo;
  fieldsError = fieldnames(outputLvlError);
  fieldsEnergy = fieldnames(outputLvlEnergy);
  fieldsHidden = fieldnames(outputLvlHidden);
  for ind = 1:length(fieldsError) 
    tableStruct.(fieldsError{ind}) = outputLvlError.(fieldsError{ind});
  end
  for ind = 2:length(fieldsEnergy) 
    tableStruct.(fieldsEnergy{ind}) = outputLvlEnergy.(fieldsEnergy{ind});
  end
  for ind = 2:length(fieldsHidden) 
    tableStruct.(fieldsHidden{ind}) = outputLvlHidden.(fieldsHidden{ind});
  end
  writetable(struct2table(tableStruct), sprintf('%s/lvlOutput.csv', dirName));

  % convergence plots
  convergenceFig = figure('visible', figVisible);
  loglog(nrDof4lvl, eta4lvl, '-o');
  hold on
  loglog(nrDof4lvl, etaVol4lvl, '-o');
  loglog(nrDof4lvl, etaJumps4lvl, '-o');
  convergenceFigLegend = ...
    [sprintf("\\eta"), sprintf("\\eta_{Vol}"), sprintf("\\eta_{Jumps}")];
  if exactSolutionKnown
    loglog(nrDof4lvl, error4lvl, '-o');
    convergenceFigLegend(end+1) = sprintf("||u-u_{CR}||_{L^2}");
    if useExactEnergy
      loglog(nrDof4lvl, diffDiscExactE4lvl, '-o');
      convergenceFigLegend(end+1) = sprintf("|E(u)-E_{NC}(u_{CR})|");
    end
  end
  if rhsGradientKnown
    if useExactEnergy
      loglog(nrDof4lvl, diffExactEGleb4lvl, '-o');
      convergenceFigLegend(end+1) = "E(u) - GLEB";
    end
    loglog(nrDof4lvl, diffGuebGleb4lvl, '-o');
    convergenceFigLegend(end+1) = "E_{NC}(J_1 u_{CR}) - GLEB";
  end
  legend(convergenceFigLegend, 'Location', 'SW');
  xlabel('nrDof');
  fName = sprintf('%s/convergence.png', dirName);
  saveas(convergenceFig, fName);

  % time and nrIteration development
  miscFig = figure('visible', figVisible);
  loglog(nrDof4lvl, time4lvl, '-o');
  hold on
  loglog(nrDof4lvl, nrIterations4lvl, '-o');
  legend([sprintf("time"), sprintf("number of iterations")], 'Location', 'SE');
  xlabel('nrDof');
  fName = sprintf('%s/miscInfo.png', dirName);
  saveas(miscFig, fName);
end
