function analyticExampleProlongation(red)
addpath(genpath(pwd));

alpha = 1; 
delta = 1;     

f=@(x)g(x,alpha,delta);  
uExact=@(x)gUexact(x,alpha,delta);  
    
[c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red);

n4s = computeN4s(n4e);
%mid4s = computeMid4s(c4n,n4s);
%
%u = f(mid4s);
u = interpolationNC(f,c4n,n4e,n4s);
du = computeGradientNC(c4n,n4e,u);
Lambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
Lambda(isinf(Lambda)) = 0;

%  Lambda = zeros(size(n4e,1),2);
%  u = zeros(size(n4s,1),1);

terminate = 0.0001;   
corr = terminate+1;
ItStep=0;
figVisible = 'off';
% set(0,'DefaultFigureVisible','off');

%% Prepare saving of results
% 'program' should be the running directory

%dirName = sprintf('../../tex/pictures/Experiment0002/quotientAbsVu/alpha_%d_beta_%d_initalRed_%d',alpha,delta,red);
dirName = sprintf('../../tex/pictures/experiment/alpha_%d_beta_%d_initalRed_%d',alpha,delta,red);
%dirName = sprintf('../../tex/pictures/experimentLongRun/alpha_%d_beta_%d_initalRed_%d',alpha,delta,red);
%dirName = sprintf('../../tex/pictures/Experiment0001/energyDiff/alpha_%d_beta_%d_initalRed_%d',alpha,delta,red);

warning('off','MATLAB:MKDIR:DirectoryExists');
mkdir(dirName);
warning('on','MATLAB:MKDIR:DirectoryExists');

% corrFig = figure();
% set(corrFig,'visible','off')
% hold on

enFig = figure('visible',figVisible);
hold on;

errExactFig = figure('visible',figVisible);
hold on;

legendEntry={};
legendEntryExactError={};
maxLengthEnergy = 0;

%% Main
while corr > terminate;
  h = 2^(-red-ItStep);
  epsStop = h^2/8;
  
  tic;
  [u,corr,corr_vec,energy_vec,errorExactVec] = ...
    tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,u,Lambda,f,alpha,...
    max(terminate,epsStop),uExact);
  time = toc; 
  
  % plot approximations
  approxFig = figure('visible',figVisible); 
  plotCR(c4n,n4e,u);
  ftitle=sprintf('approximation for red=%d, \\alpha =%d, \\beta =%d',...
  red+ItStep,alpha,delta);
  title(ftitle);
  fName = sprintf('%s/solution_red_%d.png',dirName,red+ItStep);
  saveas(approxFig,fName);
  
  approxFigAxis = figure('visible',figVisible); 
  plotAxisNC(c4n,n4e,u);
  ftitle=sprintf('approximation along axis for red=%d, \\alpha =%d, \\beta =%d',...
  red+ItStep,alpha,delta);
  title(ftitle);
  fName = sprintf('%s/solution_red_%d_axis.png',dirName,red+ItStep);
  saveas(approxFigAxis,fName);
  
  
  [c4nNew,n4eNew,n4sDbNew,n4sNbNew] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
  u = computeRefinementExtension(c4n,n4e,c4nNew,n4eNew,u);
  c4n = c4nNew;
  n4e = n4eNew;
  n4sDb = n4sDbNew;
  n4sNb = n4sNbNew;
  temp=unique(n4sDb);
  c4n(temp,:)=c4n(temp,:)./repmat(sqrt(c4n(temp,1).^2+c4n(temp,2).^2),1,2);
  du = computeGradientNC(c4n,n4e,u);
  Lambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
  Lambda(isinf(Lambda)) = 0;

  %% Plots    
  %corrFig = figure('visible',figVisible);
  %loglog(corr_vec);
  %ftitle=sprintf('loglog plot - corr for red=%d, \\alpha =%d, \\beta =%d',...
  %    red+ItStep,alpha,delta);
  %title(ftitle);
  %xlabel('number of iterations');
  %ylabel('corr');
  %fName = sprintf('%s/corr_red_%d_loglog.png',dirName,red+ItStep);
  %saveas(corrFig,fName);
  %
  %corrFig = figure('visible',figVisible);
  %plot(corr_vec);
  %ftitle=sprintf('plot - corr for red=%d, \\alpha =%d, \\beta =%d',...
  %    red+ItStep,alpha,delta);
  %title(ftitle);
  %xlabel('number of iterations');
  %ylabel('corr');
  %fName = sprintf('%s/corr_red_%d.png',dirName,red+ItStep);
  %saveas(corrFig,fName);

      
  figure(enFig);
  plot(energy_vec);
  hold on
  legendEntry(end+1)={sprintf('red = %d (%0.2fs)',red+ItStep,time)};
  if length(energy_vec) > maxLengthEnergy
    maxLengthEnergy = length(energy_vec);
  end
  
  figure(errExactFig);
  plot(errorExactVec);
  hold on
  legendEntryExactError(end+1)={sprintf('red = %d (%0.2fs)',red+ItStep,time)};

  enDiffExactFig = figure('visible',figVisible);
  loglog(abs(energy_vec+2.05802391003896));
  ftitle=sprintf('loglog plot - |E_{NC}(u_{NC})-E_u| for red=%d, \\alpha =%d, \\beta =%d',...
      red+ItStep,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{NC})-E_u|');
  fName = sprintf('%s/enDiffExact_red_%d.png',dirName,red+ItStep);
  saveas(enDiffExactFig,fName);

  singleEnergy = figure('visible',figVisible);
  plot(energy_vec);
  hold on
  plot(-2.05802391003896*ones(1,length(energy_vec)));
  ftitle=sprintf('plot - energy development for red=%d, \\alpha =%d, \\beta =%d',...
      red+ItStep,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('E_{NC}(u_{NC})');
  fName = sprintf('%s/singleEnergy_red_%d.png',dirName,red+ItStep);
  saveas(singleEnergy,fName);

  singleExactError = figure('visible',figVisible);
  loglog(errorExactVec);
  ftitle=sprintf('loglog - exact error for red=%d, \\alpha =%d, \\beta =%d',...
      red+ItStep,alpha,delta);
  title(ftitle);
  xlabel('number of iterations');
  ylabel('|E_{NC}(u_{NC})-E_u|');
  fName = sprintf('%s/singleExactError%d.png',dirName,red+ItStep);
  saveas(singleExactError,fName);

  ItStep=ItStep+1;
end 

figure(enFig);
plot(-2.05802391003896*ones(1,maxLengthEnergy));
legendEntry(end+1)={sprintf('E_u=-2.05802391003896')};
ftitle=sprintf('Energy for inital red=%d, \\alpha =%d, \\beta =%d',...
red,alpha,delta);
title(ftitle);
xlabel('number of iterations');
ylabel('energy');
legend(legendEntry);
fName = sprintf('%s/energy.png',dirName);
saveas(gcf,fName);

figure(errExactFig);
ftitle=sprintf('Exact error for inital red=%d, \\alpha =%d, \\beta =%d',...
red,alpha,delta);
title(ftitle);
xlabel('number of iterations');
ylabel('exact Error');
legend(legendEntryExactError);
fName = sprintf('%s/exactError.png',dirName);
saveas(gcf,fName);

function generateOutput(S, caseNr, k)
    name = ['../data/lshape',num2str(caseNr),'_',num2str(k),'.csv'];
    file = fopen(name,'w');
    fprintf(file, 'a,b,c,d,e,f,g,h\r\n');
    maxLvl = size(S.level,2);
    for lvl = 1:maxLvl
       fprintf(file, '%d,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%d\r\n',2^lvl,...
           S.level(lvl).errH1,S.level(lvl).rateH1,...
           S.level(lvl).errL2,S.level(lvl).rateL2,S.level(lvl).errEnergy,...
           S.level(lvl).eta,S.level(lvl).ndof);
    end
    fclose(file);
end
