function example(red,terminate)
addpath(genpath(pwd));

alpha = 1; 
delta = 1;     

%% given analytic example
f=@(x)g(x,alpha,delta);  
uExact=@(x)gUexact(x,alpha,delta);  

%% f = 0
%f=@(x)0;  
%uExact=@(x)0;  

    
[c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red);

n4s = computeN4s(n4e);
%mid4s = computeMid4s(c4n,n4s);
%
%u = f(mid4s);

u = interpolationNC(f,c4n,n4e,n4s);
du = computeGradientNC(c4n,n4e,u);
Lambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
Lambda(isinf(Lambda)) = 0;
Lambda(isnan(Lambda)) = 0;

%  Lambda = zeros(size(n4e,1),2);
%  u = zeros(size(n4s,1),1);

figVisible = 'off';
% set(0,'DefaultFigureVisible','off');

%% Prepare saving of results
% 'program' should be the running directory

dirName = sprintf('../../tex/pictures/experiment/fEqualsZero/alpha_%d_beta_%d_initalRed_%d_terminate_%f',alpha,delta,red,terminate);

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

tic;
[u,corr,corr_vec,energy_vec,errorExactVec] = ...
  tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,u,Lambda,f,alpha,...
  terminate,uExact);
time = toc; 

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
legendEntry(end+1)={sprintf('red = %d (%0.2fs)',red,time)};
if length(energy_vec) > maxLengthEnergy
  maxLengthEnergy = length(energy_vec);
end

figure(errExactFig);
plot(errorExactVec);
hold on
legendEntryExactError(end+1)={sprintf('red = %d (%0.2fs)',red,time)};

enDiffExactFig = figure('visible',figVisible);
loglog(abs(energy_vec+2.05802391003896));
ftitle=sprintf('loglog plot - |E_{NC}(u_{NC})-E_u| for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
title(ftitle);
xlabel('number of iterations');
ylabel('|E_{NC}(u_{NC})-E_u|');
fName = sprintf('%s/enDiffExact_red_%d.png',dirName,red);
saveas(enDiffExactFig,fName);

singleEnergy = figure('visible',figVisible);
plot(energy_vec);
hold on
plot(-2.05802391003896*ones(1,length(energy_vec)));
ftitle=sprintf('plot - energy development for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
title(ftitle);
xlabel('number of iterations');
ylabel('E_{NC}(u_{NC})');
fName = sprintf('%s/singleEnergy_red_%d.png',dirName,red);
saveas(singleEnergy,fName);

singleExactError = figure('visible',figVisible);
loglog(errorExactVec);
ftitle=sprintf('loglog - exact error for red=%d, \\alpha =%d, \\beta =%d',...
    red,alpha,delta);
title(ftitle);
xlabel('number of iterations');
ylabel('|E_{NC}(u_{NC})-E_u|');
fName = sprintf('%s/singleExactError%d.png',dirName,red);
saveas(singleExactError,fName);




% Plots
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
