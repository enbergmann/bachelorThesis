function exp1(red,terminate)
  
  addpath(genpath('../'),genpath('../../../utils/'));

  figVisible = 'off';

  initalU = 'zero'; %f, zero

  alpha = 1;
  delta = 1;

  h = 2^(-red); 
  tau = h^(1/2)/10; 
  
  [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(red); 

  f = @(x)g(x,alpha,delta);  
  uExact = @(x)gUexact(x,alpha,delta);

  if strcmp(initalU,'zero')
    nC = size(c4n,1);
    u = zeros(nC,1); 
    message = sprintf('on unit circle, inital u = 0');
    dirInfoName = sprintf('zeroInitial');
  elseif strcmp(initalU,'f')
    u = f(c4n);
    message = sprintf('on unit circle, inital u = I_NC(f)');
    dirInfoName = sprintf('fInitial');
  end
  
  tic;
  [corrVec,energyVec,u] = ...
    tvRegPrimalDual(c4n,n4e,n4sDb,n4sNb,h,tau,red,terminate,alpha,f,u);
  time = toc;
  
  saveResults('S1','exp1',dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
    terminate,time,corrVec,energyVec);
end
