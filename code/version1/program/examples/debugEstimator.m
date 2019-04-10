function debugEstimator

  addpath(genpath('../'),genpath('../../../utils/'));


  alpha = 1; 
  delta = 1;     
  
  [c4n,n4e,n4sDb,n4sNb] = computeGeometryCriss(0);
  n4s = computeN4s(n4e);

  f=@(x)1;  
  u = [0;1;2;3;4];

  eta4e = estimateError4e(u,f,c4n,n4e,n4sDb,n4sNb,alpha,delta)
  eta4nrDoF = sum(eta4e);
end
