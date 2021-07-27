function [c4n,n4e,u]=tv_reg_primal_dual(red,terminate)
% Input red: Anzahl der durchzuführenden Rotverfeinerungen

% Die benötigten Parameter werden bereitgestellt.
  addpath(genpath(pwd));

  alpha = 1;
  delta = 1;
  h = 2^(-red); 
  tau = h^(1/2)/10; 
  
  format long;
  % Die Geometrie wird erzeugt.
  [c4n,n4e,n4sDb,~] = computeGeometryPolygon(red); 
  
  % Die Steifigkeits- und Massematrix werden erzeugt.
  [s,m] = fe_matrices(c4n,n4e);
  
  % Die gemischte Matrix wird erzeugt.
  ms = mixed_matrix(c4n,n4e);
  
  %
  A = m+h*s; 
  
  nC = size(c4n,1); % Anzahl der Knoten
  nE = size(n4e,1); % Anzahl der Elemente
  
  % Die Funktionswerte an den Knoten werden berechnet.
  f = @(x)g(x,alpha,delta);  
  uExact = @(x)gUexact(x,alpha,delta);

  % Startwerte
  u = zeros(nC,1); 

  % u = f(c4n);
  u_tilde = u; 
  p = zeros(nE,2); 

  message = sprintf('on unit circle, inital u = 0');
  % message = sprintf('on unit circle, inital u = I_NC(f)');

  dirInfoName = sprintf('zeroInitial');
  % dirInfoName = sprintf('fInitial');
  
  %%
  nrElems = size(n4e,1);
  area4e = computeArea4e(c4n,n4e);
  s4e = computeS4e(n4e);
  nrSides = max(max(s4e));
  [~,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);
  n4s = computeN4s(n4e);

  %% Dirichlet boundary conditions
  DbNodes = unique(n4sDb);           % Dirichlet boundary nodes
  dof = setdiff(1:nC,DbNodes); % free nodes to be approximated
  % x = zeros(nC,1); % get the Dirichlet values
  % DbCoords = c4n(DbNodes,:); % coordinates of Dirichlet nodes
  % x(DbNodes) = u4Db(DbCoords);
  % b = b - A * x;  % substract inhomogenous boundary


  tic;
  [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);
  temp = zeros(nrSides,1);
  for elem = 1:1:nrElems
    temp(s4e(elem,:)) = temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
  end

  corr = 1; % ||d_t u_h^k||_{h,1/2} 
  %terminate = 0.01;
  corr_vec=[];
  energy_vec=[];
  errorExactVec=[];
  step=0;
  E = 1;
  while corr > terminate
      step=step+1;
      du_tilde = comp_gradient(c4n,n4e,u_tilde);
      p_tmp = p+tau*du_tilde;
      p = p_tmp./max(1,(sqrt(sum(p_tmp.^2,2))*ones(1,2))); 
      P = reshape(p',2*nE,1);
      
      u_new = zeros(nC,1);
      LHS = (A+tau*alpha*m);
      RHS = (A*u-tau*ms*P+tau*alpha*m*f(c4n));
      % u_new(dof) = LHS(dof,dof)\RHS(dof);
      u_new = LHS\RHS;
      dt_u = (u-u_new)/tau; 
      u_tilde = 2*u_new-u; 
      
      %corr = sqrt(dt_u'*A*dt_u);
      uNC = convertS1toCR(n4s,u);
      duNC = computeGradientNC(c4n,n4e,uNC);
      ENew = computeEnergy(area4e,uNC,duNC,alpha,temp,MAMANC);
      corr = abs(ENew-E);
      E = ENew;

      fprintf('corr/epsStop: %e / %e\n',corr,terminate);
      format long;
      fprintf('E = %f, E_exact = %f\n', E, -2.05802391003896);
      format short;
      fprintf('============================== \n');

      energy_vec=cat(2,energy_vec,ENew);
      corr_vec=cat(2,corr_vec,corr);
      errorExactVec=cat(2,errorExactVec,...
        sqrt(sum(error4eP1L2(c4n,n4e,uExact,u_new))));

      u = u_new;
      trisurf(n4e,c4n(:,1),c4n(:,2),u);
      drawnow;
  end
  time = toc;
  

  %% Prepare saving of results
  % 'program' should be the running directory
  
  dirName = sprintf('../../tex/results/conforming/%s/%s',...
    dirInfoName,datestr(now,'yy_mm_dd_HH_MM_SS'));
  
  warning('off','MATLAB:MKDIR:DirectoryExists');
  mkdir(dirName);
  warning('on','MATLAB:MKDIR:DirectoryExists');
  
  name = sprintf('%s/workspace.mat',dirName);
  save(name);

  % plot approximations
  figVisible = 'on';

  approxFig = figure('visible',figVisible); 
  trisurf(n4e,c4n(:,1),c4n(:,2),u);
  ftitle=sprintf('approximation for red=%d, \\alpha =%d, \\beta =%d, \\epsilon_{stop}=%2.3f',...
      red,alpha,delta,terminate);
  title(ftitle);
  fName = sprintf('%s/solution_red_%d.png',dirName,red);
  saveas(approxFig,fName);
  
  approxFigAxis = figure('visible',figVisible); 
  plotAxis(c4n,u);
  ftitle=sprintf('approximation along axis for red=%d, \\alpha =%d, \\beta =%d, \\epsilon_{stop}=%2.3f',...
      red,alpha,delta,terminate);
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
  fprintf(file, '%.8e\n',corr_vec);
  fclose(file);
  
  name = sprintf('%s/energyVec.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, '%.8g\n',energy_vec);
  fclose(file);
  
  name = sprintf('%s/errorExactVec.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, '%.8e\n',errorExactVec);
  fclose(file);
  
  name = sprintf('%s/exactAbsoluteEnergyDifference.txt',dirName);
  file = fopen(name,'w');
  fprintf(file, '%.8g\n',abs(energy_vec+2.05802391003896));
  fclose(file);
  
  % further plots  
  enFig = figure('visible',figVisible); 
  loglog(energy_vec);
  hold on;
  plot(-2.05802391003896*ones(1,length(energy_vec)));
  legend(sprintf('red = %d (%0.2fs)',red,time),sprintf('E_u=-2.05802391003896'));
  hold off;
  ftitle=sprintf('Energy for inital red=%d, \\alpha =%d, \\beta =%d',...
  red,alpha,delta);
  title(ftitle);
  fName = sprintf('%s/energy.png',dirName);
  xlabel('number of iterations');
  ylabel('energy');
  saveas(enFig,fName);
  
  errExactFig = figure('visible',figVisible);
  loglog(errorExactVec);
  ftitle=sprintf('Exact error for inital red=%d, \\alpha =%d, \\beta =%d',...
  red,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('exact Error');
  legend(sprintf('red = %d (%0.2fs)',red,time));
  fName = sprintf('%s/exactError.png',dirName);
  saveas(gcf,fName);
  
  enDiffExactFig = figure('visible',figVisible);
  loglog(abs(energy_vec+2.05802391003896));
  ftitle=sprintf('|E_{NC}(u_{NC})-E_u| for red=%d, \\alpha =%d, \\beta =%d',...
      red,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{NC})-E_u|');
  fName = sprintf('%s/enDiffExact_red_%d.png',dirName,red);
  saveas(enDiffExactFig,fName);
  
  corrFig = figure('visible',figVisible);
  loglog(corr_vec);
  ftitle=sprintf('loglog plot - corr for red=%d, \\alpha =%d, \\beta =%d',...
      red,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('corr');
  fName = sprintf('%s/corr_red_%d_loglog.png',dirName,red);
  saveas(corrFig,fName);
end
