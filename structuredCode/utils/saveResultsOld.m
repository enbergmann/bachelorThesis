%TODO probably folder documentation or sth for saveRestulst, saveScreenshots and
%so on

%write two of these, one for CR and one for S1

%function saveResults(fe,experiment,dirInfoName,figVisible, ...
%    message,c4n,n4e,u,red,alpha,delta, ...
%    terminate,time,corrVec,energyVec,tau,miscMsg,nrDof,...
%    useAdaptivity,nrDoF4lvl,eta4lvl,l2Error4lvl)

function saveResults(fe,experiment,dirInfoName,figVisible, ...
    message,c4n,n4e,u,red,alpha,delta, ...
    terminate,time,corrVec,energyVec,tau,miscMsg,nrDof,...
    useAdaptivity,nrDoF4lvl,eta4lvl,l2Error4lvl)

  warning('off','MATLAB:MKDIR:DirectoryExists');
  if strcmp(fe,'S1')
    dirName = sprintf('../../../../results/%s/conforming/%s/%s',...
      experiment,dirInfoName,datestr(now,'yy_mm_dd_HH_MM_SS'));
    mkdir(dirName);
    approxFig = figure('visible',figVisible); 
    show_p1(c4n,n4e,u);
    ftitle=sprintf('approximation for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
    title(ftitle);
    fName = sprintf('%s/solution_red_%d.png',dirName,red);
    saveas(approxFig,fName);
    
    approxFigAxis = figure('visible',figVisible); 
    plotAxis(c4n,u);
    ftitle=sprintf('approximation along axis for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
    title(ftitle);
    fName = sprintf('%s/solution_red_%d_axis.png',dirName,red);
    saveas(approxFigAxis,fName);
  elseif strcmp(fe,'CR')
    dirName = sprintf('../../../../results/%s/nonconforming/%s/%s',...
      experiment,dirInfoName,datestr(now,'yy_mm_dd_HH_MM_SS'));
    mkdir(dirName);
    approxFig = figure('visible',figVisible); 
%    plotCR(c4n,n4e,u);
    plotCR(c4n,n4e,u,{'CR Solution'; [num2str(nrDof) ' degrees of freedom']});
    ftitle=sprintf('approximation for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
    title(ftitle);
    fName = sprintf('%s/solution_red_%d.png',dirName,red);
    saveas(approxFig,fName);
    
    approxFigAxis = figure('visible',figVisible); 
    plotAxisNC(c4n,n4e,u);
    ftitle=sprintf('approximation along axis for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
    title(ftitle);
    fName = sprintf('%s/solution_red_%d_axis.png',dirName,red);
    saveas(approxFigAxis,fName);
  else
    dirName = sprintf('../../../results/unknownFE/%s/%s',...
      dirInfoName,datestr(now,'yy_mm_dd_HH_MM_SS'));
    feKnown = 0;
    mkdir(dirName);
  end
  warning('on','MATLAB:MKDIR:DirectoryExists');
    
  name = sprintf('%s/workspace.mat',dirName);
  save(name,'-regexp','^(?!(c4n|n4e|approxFig|approxFigAxis)$).');
  % save all except for the listed in the end
   
  name = sprintf('%s/setting.txt',dirName);
  file = fopen(name,'w');
  fprintf(file,...
    '%s\nalpha = %.8g \nbeta = %.8g \nred = %d \nnrDof = %d \nepsStop = %.2e\ntime = %0.2fs\ntau = %.8g\n\nmisc: %s\n'...
    ,message,alpha, delta, red, nrDof, terminate, time, tau, miscMsg);
  fclose(file);
  
  name = sprintf('%s/corrVec.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, '%.8e\n',corrVec);
  fclose(file);
  
  name = sprintf('%s/energyVec.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, '%.8g\n',energyVec);
  fclose(file);
  
  name = sprintf('%s/exactAbsoluteEnergyDifference.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, '%.8g\n',abs(energyVec+2.05802391003896));
  fclose(file);
  
  % further plots  
  enFig = figure('visible',figVisible); 
  plot(energyVec);
  hold on;
  plot(-2.05802391003896*ones(1,length(energyVec)));
  legend(sprintf('red = %d (%0.2fs)',red,time),sprintf('E_u=-2.05802391003896'));
  hold off;
  ftitle=sprintf('Energy for inital red=%d, \\alpha =%d, \\beta =%d',...
  red,alpha,delta);
  title(ftitle);
  fName = sprintf('%s/energy.png',dirName);
  xlabel('number of iterations');
  ylabel('energy');
  saveas(enFig,fName);
  
  enDiffExactFig = figure('visible',figVisible);
  loglog(abs(energyVec+2.05802391003896));
  ftitle=sprintf('|E_{NC}(u_{NC})-E_u| for red=%d, \\alpha =%d, \\beta =%d',...
      red,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{NC})-E_u|');
  fName = sprintf('%s/enDiffExact_red_%d.png',dirName,red);
  saveas(enDiffExactFig,fName);
  
  corrFig = figure('visible',figVisible);
  loglog(corrVec);
  ftitle=sprintf('loglog plot - corr for red=%d, \\alpha =%d, \\beta =%d',...
      red,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('corr');
  fName = sprintf('%s/corr_red_%d_loglog.png',dirName,red);
  saveas(corrFig,fName);

  triangFig = figure('visible',figVisible);
  plotTriangulation(c4n,n4e);
  ftitle = sprintf('Triangulation for red=%d, \\alpha =%d, \\beta =%d',...
      red,alpha,delta);
  title(ftitle);
  fName = sprintf('%s/triangulation.png',dirName);
  saveas(triangFig,fName);

  if useAdaptivity
    estimatorAndL2ErrorFig = figure('visible',figVisible);
    % plotConvergence(nrDoF4lvl,eta4lvl,'\eta_l');
    loglog(nrDoF4lvl,eta4lvl);
    hold on
    loglog(nrDoF4lvl,l2Error4lvl);
    legend(sprintf('\\eta'),sprintf('||u-u_{CR}||_{L^2}'));
    ftitle = sprintf('Error and Estimator for red=%d, \\alpha =%d, \\beta =%d',...
        red,alpha,delta);
    title(ftitle);
    xlabel('number of iterations');
    ylabel('corr');
    fName = sprintf('%s/estimator.png',dirName);
    saveas(estimatorAndL2ErrorFig,fName);

    name = sprintf('%s/nrDoF4lvl.txt',dirName);
    file = fopen(name,'w');
    fprintf(file, '%.8g\n',nrDoF4lvl);
    fclose(file);

    name = sprintf('%s/eta4lvl.txt',dirName);
    file = fopen(name,'w');
    fprintf(file, '%.8g\n',eta4lvl);
    fclose(file);

    name = sprintf('%s/l2Error4lvl.txt',dirName);
    file = fopen(name,'w');
    fprintf(file, '%.8g\n',l2Error4lvl);
    fclose(file);
  end
end
