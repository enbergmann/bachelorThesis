% TODO probably folder called documentation or sth for saveResults,
%      saveScreenshots and so on
% TODO this is all not done yet, just changed so it functions again

%function saveResults(experiment,dirInfoName,figVisible, ...
%    message,c4n,n4e,u,red,alpha,delta, ...
%    terminate,time,corrVec,energyVec,tau,miscMsg,nrDof,...
%    useAdaptivity,nrDof4lvl,eta4lvl,error4lvl)

function saveResults(params, currData, output)

  % extract necessary parameters from params
  expName = params.expName;
  dirInfoName = params.dirInfoName;
  figVisible = params.figVisible;
  parAlpha = params.parAlpha;
  parBeta = params.parBeta;

  % extract necessary information from currData
  c4n = currData.c4n;
  n4e = currData.n4e;
  nrDof = currData.nrDof; 

  % extract necessary information from output
  u = output.u;
  nrDof4lvl = output.nrDof4lvl;
  corrVec = output.corrVec;
  energyVec = output.energyVec;
  time = output.time;
  nrDof4lvl = output.nrDof4lvl;
  eta4lvl = output.eta4lvl;
  error4lvl = output.error4lvl;

  % create directory
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  % startAlgorithmNC run from .../structuredcode/nonconforming
  dirName = sprintf('../../results/nonconforming/%s/%s/%s', ...
    expName, dirInfoName, datestr(now,'yy_mm_dd_HH_MM_SS'));
  mkdir(dirName);
  warning('on', 'MATLAB:MKDIR:DirectoryExists');
    
  % TODO think about what stuff might be relevant from the given structs and
  %      save them if there is anything useful
  % % save workspace
  % name = sprintf('%s/workspace.mat', dirName);
  % save(name,'-regexp','^(?!(c4n|n4e|approxFig|approxFigAxis)$).');
  % % save all except for the listed in the end
  
  % save results
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

  % TODO think about what should be saved in setting.txt, maybe one can 
  %      easily save everything from struct (google it)
  % name = sprintf('%s/setting.txt', dirName);
  % file = fopen(name, 'w');
  % fprintf(file, ...
  %   ['%s\nalpha = %.8g \nbeta = %.8g \nnrDof = %d ', 
  %   '\nepsStop = %.2e\ntime = %0.2fs\ntau = %.8g\n\nmisc: %s\n'], ...
  %   message, parAlpha, parBeta, nrDof, terminate, time, tau, miscMsg);
  % fclose(file);
  
  name = sprintf('%s/corrVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8e\n', corrVec);
  fclose(file);
  
  name = sprintf('%s/energyVec.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', energyVec);
  fclose(file);
  
  % TODO THIS IS ALL JUST A TEST, REWRITE AND CHANGE DETAILS AND INFOS
  %TODO here the given (!) exact energy (if exactEnergyKnow) should be used
  %     also from here on out every time it's used
  %     this means it still needs to be done how the externly calculated 
  %     energy should be communicated to the code (prob see editable)
  name = sprintf('%s/exactAbsoluteEnergyDifference.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', abs(energyVec+2.05802391003896));
  fclose(file);
  
  % further plots  
  enFig = figure('visible', figVisible); 
  plot(energyVec);
  hold on;
  plot(-2.05802391003896*ones(1, length(energyVec)));
  legend(sprintf('nrDof = %d (%0.2fs)', nrDof, time), ...
    sprintf('E_u=-2.05802391003896'));
  hold off;
  ftitle=sprintf('Energy for inital nrDof=%d, \\alpha =%d, \\beta =%d',...
  nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/energy.png', dirName);
  xlabel('number of iterations');
  ylabel('energy');
  saveas(enFig, fName);
  
  enDiffExactFig = figure('visible',figVisible);
  loglog(abs(energyVec+2.05802391003896));
  ftitle = sprintf('|E_{NC}(u_{NC})-E_u| for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{NC})-E_u|');
  fName = sprintf('%s/enDiffExact_red_%d.png', dirName, nrDof);
  saveas(enDiffExactFig, fName);
  
  corrFig = figure('visible', figVisible);
  loglog(corrVec);
  ftitle=sprintf('loglog plot - corr for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('corr');
  fName = sprintf('%s/corr_nrDof_%d_loglog.png', dirName, nrDof);
  saveas(corrFig, fName);

  triangFig = figure('visible',figVisible);
  plotTriangulation(c4n,n4e);
  ftitle = sprintf('Triangulation for nrDof = %d, \\alpha = %d, \\beta = %d', ...
    nrDof, parAlpha, parBeta);
  title(ftitle);
  fName = sprintf('%s/triangulation.png', dirName);
  saveas(triangFig, fName);

  % % TODO careful, what error is the estimator for (this error should always
  % % be computed then, make error4lvl a alternative then)
  % % TODO use exact solutionKnown if want to use error4lvl, else it's empty
  % estimatorAndErrorFig = figure('visible', figVisible);
  % % plotConvergence(nrDof4lvl,eta4lvl,'\eta_l');
  % loglog(nrDof4lvl, eta4lvl);
  % hold on
  % loglog(nrDof4lvl, error4lvl);
  % legend(sprintf('\\eta'), sprintf('||u-u_{CR}||_{L^2}'));
  % ftitle = sprintf('Error and Estimator for nrDof = %d, \\alpha = %d, \\beta = %d', ...
  %     nrDof, parAlpha, parBeta);
  % title(ftitle);
  % xlabel('number of iterations');
  % ylabel('corr');
  % fName = sprintf('%s/estimator.png', dirName);
  % saveas(estimatorAndErrorFig, fName);

  name = sprintf('%s/nrDof4lvl.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', nrDof4lvl);
  fclose(file);

  name = sprintf('%s/eta4lvl.txt', dirName);
  file = fopen(name, 'w');
  fprintf(file, '%.8g\n', eta4lvl);
  fclose(file);

  % name = sprintf('%s/error4lvl.txt', dirName);
  % file = fopen(name, 'w');
  % fprintf(file, '%.8g\n', error4lvl);
  % fclose(file);
end
