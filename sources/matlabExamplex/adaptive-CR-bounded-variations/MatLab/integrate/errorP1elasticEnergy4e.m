function [elastErrorSq4e,elastEnergyUexSq4e,elastEnergyUsq4e] = ...
    errorP1elasticEnergy4e(geometryData,exSolutionData,localData,solutionData,parameter)
  lame_mu = parameter.lame_mu;
  lambda = parameter.lambda;
  intDegree = parameter.intDegree;
  % works for any piecewise P1 function (CR and KS) included
  %% LOAD DATA FROM CELL-ARRAYs
  c4n     = geometryData.c4n; %geometryData{1};
  n4e     = geometryData.n4e; %geometryData{2};
  n4sDb1C = geometryData.n4sDb1C; %geometryData{3};
  n4sNb1C = geometryData.n4sNb1C; %geometryData{4};
  n4sDb2C = geometryData.n4sDb2C; %geometryData{5};
  n4sNb2C = geometryData.n4sNb2C; %geometryData{6};
  
  area4e  = localData.area4e; %localData{15};
  
  sigmaExact = exSolutionData.sigmaExact; %exSolutionData{7};
  epsUex  = exSolutionData.epsUex; %exSolutionData{11};
  
  epsU4e  = solutionData.epsU4e; %solutionData{9};
  
  %% ELASTICENERGY OF EXACT SOLUTION
  nrElems = size(n4e,1);
  elastEnergyUexSq4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                         permute((sum(sum(epsUex(Gpts4p).*sigmaExact(Gpts4p)))),[3,1,2]),...
                         intDegree+1,[nrElems,1,1],area4e);
  
  %% MIXED TERM
  sigmaExTmp = @(x) permute(sigmaExact(x),[3,1,2]);
  intSigma4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                            sigmaExTmp(Gpts4p),...
                            intDegree+1,[nrElems,2,2],area4e);
  intSigma4e = permute(intSigma4e,[2,3,1]);
  mixedTerm4e = sum(sum(epsU4e.*intSigma4e));
  mixedTerm4e = permute(mixedTerm4e,[3,2,1]);
  
  %% ELASTICENERGY OF P1-FUNCTION
  elastEnergyUsq4e = computeElastEnergyP1(geometryData,localData,solutionData,...
                                          lame_mu,lambda,intDegree);
  %% ADD TOGETHER
  elastErrorSq4e = elastEnergyUexSq4e - 2*mixedTerm4e + elastEnergyUsq4e;
end
