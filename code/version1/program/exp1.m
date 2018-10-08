function exp1(red,terminate)

  addpath(genpath(pwd),'../../utils/');

  figVisible = 'off';
  % set(0,'DefaultFigureVisible','off');
  
  % initalU = 'zero'; 
  initalU = 'f'; 
  
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
      

  
  %%
  
  %% Main
  
  tau = 1/2;
  h = 2^(-red);

  tic;
  [u,corrVec,energyVec] = ...
    tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,terminate,alpha,f,u,Lambda);
  time = toc; 
  
  saveResults('CR',dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
    terminate,time,corrVec,energyVec);
end
