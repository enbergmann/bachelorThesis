function exp1(red,terminate)

  addpath(genpath(pwd),'../../utils/');
  
  alpha = 1; 
  delta = 1;     
  
  [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red);
  
  n4s = computeN4s(n4e);
  
  %% given analytic example

  f=@(x)g(x,alpha,delta);  
  % uExact=@(x)gUexact(x,alpha,delta);  
  
  %% f = 0
  %f=@(x)0;  
  %uExact=@(x)0;  
  
  %mid4s = computeMid4s(c4n,n4s);
  %

  % initalU = 'zero'; 
  initalU = 'f'; 

  if strcmp(initalU,'zero')
    u = zeros(size(n4s,1),1);
    message = sprintf('on unit circle, inital u = 0');
    dirInfoName = sprintf('zeroInitial');
  elseif strcmp(initalU,'f')
    u = interpolationNC(f,c4n,n4e,n4s);
    message = sprintf('on unit circle, inital u = I_NC(f)');
    dirInfoName = sprintf('fInitial');
  end
  
  du = computeGradientNC(c4n,n4e,u);
  Lambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
  Lambda(isinf(Lambda)) = 0;
  Lambda(isnan(Lambda)) = 0;
  
  %  Lambda = zeros(size(n4e,1),2);
  %  u = zeros(size(n4s,1),1);
      

  figVisible = 'off';
  % set(0,'DefaultFigureVisible','off');
  
  %%
  
  %% Main
  
  tau = 1/2;
  h = 2^(-red);

  tic;
  [u,corrVec,energyVec] = ...
    tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,terminate,alpha,f,u,Lambda);
  time = toc; 
  
  saveResults();
  %% Prepare saving of results
  
  dirName = sprintf('../../../results/nonconforming/red%d',...
    red);
  
  warning('off','MATLAB:MKDIR:DirectoryExists');
  mkdir(dirName);
  warning('on','MATLAB:MKDIR:DirectoryExists');
  
  name = sprintf('%s/workspace.mat',dirName);
  save(name);
   
  % plot approximations
  approxFig = figure('visible',figVisible); 
  plotCR(c4n,n4e,u);
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
  
  name = sprintf('%s/setting.txt',dirName);
  file = fopen(name,'w');
  fprintf(file,...
    '%s\nalpha = %.8g \nbeta = %.8g \nred = %d \nepsStop = %.2e\ntime = %0.2fs\n'...
    ,message,alpha, delta, red, terminate, time);
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
end
