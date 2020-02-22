% TODO this is all not done yet, just changed so it functions again
%
% TODO similar to struct 2 table write those information in a file so
% one can see everything in one go, also elapsed time of program and so on

function saveResults(params, currData, outputLvlInfo, outputLvl, output)

%% INITIALIZATION
  % extract necessary parameters from params
  benchmark = params.benchmark;
  exactSolutionKnown = params.exactSolutionKnown;
  expName = params.expName;
  dirInfoName = params.dirInfoName;
  figVisible = params.figVisible;
  parAlpha = params.parAlpha;
  parBeta = params.parBeta;
  useExactEnergy = params.useExactEnergy;
  exactEnergy = params.exactEnergy;

  % extract necessary information from currData
  nrDof = currData.nrDof; 
  c4n = currData.c4n;
  n4e = currData.n4e;
  s4e = currData.s4e;

  % extract necessary information from outputLvlInfo
  currLvl = outputLvlInfo.lvl(end);
  nrDof4lvl = outputLvlInfo.nrDof;
  time = outputLvlInfo.time(end);

  % extract necessary information from outputLvlInfo
  eta4lvl = outputLvl.eta;
  etaVol4lvl = outputLvl.etaVol;
  etaJumps4lvl = outputLvl.etaJumps;
  if exactSolutionKnown
    error4lvl = outputLvl.error4lvl;
  end
  if useExactEnergy
    gleb4lvl = outputLvl.gleb;
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
  if currLvl == 0
  % this means the benchmark-file should not be changed until level 0 is saved
    source = sprintf('benchmarks/%s.m', benchmark);
    destination = ...
      sprintf('../../results/nonconforming/%s/%s/benchmark_%s.txt', ...
      expName, dirInfoName, benchmark);
    copyfile(source, destination);
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
  plotCR(c4n, n4e, u, {'CR Solution'; [num2str(nrDof) ' degrees of freedom']});
  ftitle = sprintf(...
    'approximation for nrDof = %d, \\alpha =%d, \\beta =%d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/solution_nrDof_%d.png', dirName, nrDof);
  saveas(approxFig, fName);
  
  approxFigAxis = figure('visible', figVisible); 
  plotAxisNC(c4n,n4e,u);
  ftitle = sprintf(...
    'approximation along axis for nrDof = %d, \\alpha =%d, \\beta =%d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/solution_nrDof_%d_axis.png', dirName, nrDof);
  saveas(approxFigAxis, fName);

  grayscaleFig = figure('visible', figVisible); 
  plotGrayscale(c4n, n4e, mean(u(s4e), 2), ...
    {'CR Solution'; [num2str(nrDof) ' degrees of freedom']});
  ftitle = sprintf(...
    'grayscale of approximation for nrDof = %d, \\alpha =%d, \\beta =%d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/grayscale_nrDof_%d.png', dirName, nrDof);
  saveas(grayscaleFig, fName);

%% SAVE PLOTS AND RESULTS OF THE ITERATION FOR THE LEVEL
  name = sprintf('%s/corrVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', corrVec);
  fclose(file);
  
  name = sprintf('%s/energyVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', energyVec);
  fclose(file);
  
  % TODO THIS IS ALL JUST A TEST, REWRITE AND CHANGE DETAILS AND INFOS
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

    name = sprintf('%s/gleb.txt', dirName);
    file = fopen(name, 'w');
    fprintf(file, '%.8g\n', gleb4lvl);
    fclose(file);

    glebFig = figure('visible', figVisible);
    loglog(nrDof4lvl, gleb4lvl, '-o');
    legend(sprintf('GLEB'));
    ftitle = sprintf('GLEB for nrDof = %d, \\alpha = %d, \\beta = %d', ...
      nrDof, parAlpha, parBeta);
    xlabel('nrDof');
    ylabel('GLEB');
    fName = sprintf('%s/gleb.png', dirName);
    title(ftitle);
    saveas(glebFig, fName);

    name = sprintf('%s/glebExactEnergyDifference.txt', dirName);
    file = fopen(name, 'w');
    fprintf(file, '%.8g\n', exactEnergy-gleb4lvl);
    fclose(file);

    glebExactEnergyDiffereneFig = figure('visible',figVisible);
    loglog(nrDof4lvl, exactEnergy-gleb4lvl, '-o');
    ftitle = sprintf(...
      'E_u-GLEB for nrDof = %d, \\alpha = %d, \\beta = %d', ...
      nrDof, parAlpha, parBeta);
    title(ftitle);
    xlabel('nrDof');
    ylabel('E_u-GLEB');
    fName = sprintf('%s/gleb_nrDof_%d.png', dirName, nrDof);
    saveas(glebExactEnergyDiffereneFig, fName);
  end
    
  enFig = figure('visible', figVisible); 
  plot(energyVec);
  if useExactEnergy
    hold on;
    plot(exactEnergy*ones(1, length(energyVec)));
    legend(sprintf('nrDof = %d (%0.2fs)', nrDof, time), ...
      sprintf('E_u = %.8g', exactEnergy));
    hold off;
  else
    legend(sprintf('nrDof = %d (%0.2fs)', nrDof, time));
    % TODO else unnecessary, isnt there some legend Add entry
  end
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
  % TODO save c4n and n4e to tikz triangulation (see tiens example)
  triangFig = figure('visible',figVisible);
  plotTriangulation(c4n,n4e);
  ftitle = sprintf('Triangulation for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/triangulation.png', dirName);
  saveas(triangFig, fName);

  % % TODO careful, what error is the estimator for (this error should always
  % % be computed then (if possible), make error4lvl a alternative then)

  % convergence plots
  % TODO mu and xi
  estimatorAndErrorFig = figure('visible', figVisible);
  % plotConvergence(nrDof4lvl,eta4lvl,'\eta_l');
  loglog(nrDof4lvl, eta4lvl, '-bo');
  hold on
  loglog(nrDof4lvl, etaVol4lvl, '-ms');
  hold on
  loglog(nrDof4lvl, etaJumps4lvl, '-rd');
  if exactSolutionKnown
    hold on
    loglog(nrDof4lvl, error4lvl, '-g*');
    legend(sprintf('\\eta'), ...
      sprintf('\\eta_{Vol}'), sprintf('\\eta_{Jumps}'), ...
      sprintf('||u-u_{CR}||_{L^2}'));
    ftitle = sprintf('Error and Estimator for nrDof = %d, \\alpha = %d, \\beta = %d', ...
      nrDof, parAlpha, parBeta);
    ylabel('estimator and error');
    fName = sprintf('%s/errorAndEstimator.png', dirName);
    hold off;
  else
    legend(sprintf('\\eta'), ...
      sprintf('\\eta_{Vol}'), sprintf('\\eta_{Jumps}'))
    ftitle = sprintf('Estimator for nrDof = %d, \\alpha = %d, \\beta = %d', ...
      nrDof, parAlpha, parBeta);
    ylabel('estimator');
    fName = sprintf('%s/estimator.png', dirName);
  end
  title(ftitle);
  xlabel('nrDof');
  saveas(estimatorAndErrorFig, fName);

  name = sprintf('%s/nrDof4lvl.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', nrDof4lvl);
  fclose(file);

  name = sprintf('%s/eta4lvl.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', eta4lvl);
  fclose(file);

  name = sprintf('%s/etaVol4lvl.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', etaVol4lvl);
  fclose(file);

  name = sprintf('%s/etaJumps4lvl.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', etaJumps4lvl);
  fclose(file);

  if exactSolutionKnown
    name = sprintf('%s/error4lvl.txt', dirName);
    file = fopen(name, 'w');
    fprintf(file, '%.8g\n', error4lvl);
    fclose(file);
  end
end
