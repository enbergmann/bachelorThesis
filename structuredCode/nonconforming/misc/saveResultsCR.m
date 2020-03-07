% NOTE figures that might be included in latex later without title and caption
% but everything else should have this for faster overview

% TODO similar to struct 2 table write those information in a file so
% one can see everything in one go, also elapsed time of program and so on
% TODO general file with most intersting infos about the level
% one file for AFEM output for all levels (outputLvl) (add GLEB difference)
% one file for all iteration outputs (maybe, compare to old plots)
% 
% I.E. pretty much anything can be seen in benchmark.m but make a file with 
% just the relevant information with one look, like
% alpha = ...
% beta = ...
% n = ...
% ...

function saveResultsCR(params, currData, ...
    outputLvlInfo, outputLvlError, outputLvlEnergy, output)

%% INITIALIZATION
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
  useExactEnergy = params.useExactEnergy;
  exactEnergy = params.exactEnergy;
  plotGivenFunctions = params.plotGivenFunctions;
  refinementLevel4Plots = params.refinementLevel4Plots;

  % extract necessary information from currData
  nrDof = currData.nrDof; 
  c4n = currData.c4n;
  n4e = currData.n4e;
  s4e = currData.s4e;

  % extract necessary information from outputLvlInfo
  currLvl = outputLvlInfo.lvl(end);
  nrDof4lvl = outputLvlInfo.nrDof;
  time = outputLvlInfo.time(end);

  % extract necessary information from outputLvlError
  eta4lvl = outputLvlError.eta;
  etaVol4lvl = outputLvlError.etaVol;
  etaJumps4lvl = outputLvlError.etaJumps;
  if exactSolutionKnown
    error4lvl = outputLvlError.error4lvl;
  end

  % extract necessary information from outputLvlInfoEnergy
  if useExactEnergy
    gleb4lvl = outputLvlEnergy.gleb;
  end
  
  % extract necessary information from output
  u = output.u;
  corrVec = output.corrVec;
  energyVec = output.energyVec;

%% CREATE DIRECTORY
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf('../../results/nonconforming/%s/%s/lvl%d_nrDof%d', ...
    expName, dirInfoName, currLvl, nrDof);
    % startAlgorithmNC run from .../structuredcode/nonconforming
  mkdir(dirName);
  warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
%% SAVE INFORMATION ABOUT EXPERIMENT
  tableStruct = outputLvlInfo;
  fieldsError = fieldnames(outputLvlError);
  fieldsEnergy = fieldnames(outputLvlEnergy);
  for ind = 2:length(fieldsError) 
    tableStruct.(fieldsError{ind}) = outputLvlError.(fieldsError{ind});
  end
  for ind = 2:length(fieldsEnergy) 
    tableStruct.(fieldsEnergy{ind}) = outputLvlEnergy.(fieldsEnergy{ind});
  end
  writetable(struct2table(tableStruct), sprintf('%s/lvlOutput.csv', dirName));

  if currLvl == 0
    % this means the benchmark-file should not be changed until level 0 is
    % saved
    source = sprintf('benchmarks/%s.m', benchmark);
    destination = ...
      sprintf('../../results/nonconforming/%s/%s/benchmark_%s.m', ...
      expName, dirInfoName, benchmark);
    copyfile(source, destination);

    if plotGivenFunctions
      % plot rhs and grayscale image of rhs
      % TODO make polygon mode right (there is a function for this, i.e.
      % don't use load geometry)
      if strcmp(geometry, 'Polygon'), geometry = 'BigSquare'; end
      [c4nRhs, n4eRhs] = loadGeometry(geometry, refinementLevel4Plots);

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
      axis off;
      axis equal;
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
        axis off;
        axis equal;
        colormap gray;
        fName = sprintf(...
          '../../results/nonconforming/%s/%s/exactSolutionGrayscale.png', ...
          expName, dirInfoName);
        saveas(uExactGrayscaleFig, fName);
      end
    end
  end

  % TODO think about what stuff might be relevant from the given structs and
  %      save them if there is anything useful
  %
  % % save workspace
  % name = sprintf('%s/workspace.mat', dirName);
  % save(name,'-regexp','^(?!(c4n|n4e|approxFig|approxFigAxis)$).');
  % % save all except for the listed in the end
  %
  % TODO fix it maybe, what does one need here
  % name = sprintf('%s/setting.txt', dirName);
  % file = fopen(name, 'w');
  % fprintf(file, ...
  %   ['%s\nalpha = %.8g \nbeta = %.8g \nnrDof = %d ', 
  %   '\nepsStop = %.2e\ntime = %0.2fs\ntau = %.8g\n\nmisc: %s\n'], ...
  %   message, parAlpha, parBeta, nrDof, terminate, time, tau, miscMsg);
  % fclose(file);
  %
  % TODO doesn't work, can it even?
  % name = sprintf('%s/parameters.txt', dirName);
  % writetable(struct2table(params, 'AsArray', true), name);


%% SAVE PLOTS OF SOLUTION
  approxFig = figure('visible', figVisible); 
  plotCR(c4n, n4e, u);
  fName = sprintf('%s/solution_nrDof_%d.png', dirName, nrDof);
  saveas(approxFig, fName);
  
  approxFigAxis = figure('visible', figVisible); 
  plotAxisNC(c4n,n4e,u);
  fName = sprintf('%s/solution_nrDof_%d_axis.png', dirName, nrDof);
  saveas(approxFigAxis, fName);

  grayscaleFig = figure('visible', figVisible); 
  plotGrayscale(c4n, n4e, mean(u(s4e), 2));
  fName = sprintf('%s/grayscale_nrDof_%d.png', dirName, nrDof);
  saveas(grayscaleFig, fName);

%% SAVE PLOTS AND RESULTS OF THE ITERATION FOR THE LEVEL
  % TODO this just all in one file, maybe with iterationNumber as first row?
  % --> see latex tikz first and decide after
  name = sprintf('%s/corrVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', corrVec);
  fclose(file);
  
  name = sprintf('%s/energyVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', energyVec);
  fclose(file);
  
  if useExactEnergy
    name = sprintf('%s/exactAbsoluteEnergyDifference.txt', dirName);
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
    fName = sprintf('%s/enDiffExact_nrDof_%d.png', dirName, nrDof);
    saveas(enDiffExactFig, fName);

    glebFig = figure('visible', figVisible);
    loglog(nrDof4lvl, gleb4lvl, '-o');
    ftitle = sprintf('GLEB for nrDof = %d, \\alpha = %d, \\beta = %d', ...
      nrDof, parAlpha, parBeta);
    xlabel('nrDof');
    ylabel('GLEB');
    fName = sprintf('%s/gleb.png', dirName);
    title(ftitle);
    saveas(glebFig, fName);
  end
    
  enFig = figure('visible', figVisible); 
  enFigLegend = sprintf("nrDof = %d (%0.2fs)", nrDof, time);
  plot(energyVec);
  if useExactEnergy
    hold on;
    plot(exactEnergy*ones(1, length(energyVec)));
    enFigLegend(end+1) = sprintf("E_u = %.8g", exactEnergy);
  end
  legend(enFigLegend);
  ftitle=sprintf('Energy for inital nrDof=%d, \\alpha =%d, \\beta =%d',...
  nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/energy.png', dirName);
  xlabel('number of iterations');
  ylabel('energy');
  saveas(enFig, fName);
  
  
  corrFig = figure('visible', figVisible);
  loglog(corrVec);
  ftitle=sprintf('loglog plot - corr for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('corr');
  fName = sprintf('%s/corr_nrDof_%d_loglog.png', dirName, nrDof);
  saveas(corrFig, fName);

%% SAVE AFEM RESULTS AND TRIANGULATION
  triangFig = figure('visible',figVisible);
  plotTriangulation(c4n,n4e);
  fName = sprintf('%s/triangulation.png', dirName);
  saveas(triangFig, fName);

  dlmwrite(sprintf('%s/c4n.txt', dirName), c4n, 'Delimiter', '\t');
  dlmwrite(sprintf('%s/n4e.txt', dirName), n4e, 'Delimiter', '\t');

  % convergence plots 

  % % TODO careful, what error is the estimator for (this error should always
  % % be computed then (if possible), make error4lvl a alternative then)
  convergenceFig = figure('visible', figVisible);
  loglog(nrDof4lvl, eta4lvl, '-o');
  hold on
  loglog(nrDof4lvl, etaVol4lvl, '-o');
  hold on
  loglog(nrDof4lvl, etaJumps4lvl, '-o');
  convergenceFigLegend = ...
    [sprintf("\\eta"), sprintf("\\eta_{Vol}"), sprintf("\\eta_{Jumps}")];
  if exactSolutionKnown
    hold on
    loglog(nrDof4lvl, error4lvl, '-o');
    convergenceFigLegend(end+1) = sprintf("||u-u_{CR}||_{L^2}");
    if useExactEnergy
      hold on
      loglog(nrDof4lvl, exactEnergy-gleb4lvl, '-o');
      convergenceFigLegend(end+1) = "GLEB";
    end
  end
  legend(convergenceFigLegend, 'Location', 'SW');
  xlabel('nrDof');
  fName = sprintf('%s/convergence.png', dirName);
  saveas(convergenceFig, fName);
end
