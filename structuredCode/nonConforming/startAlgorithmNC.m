function [params, output] = startAlgorithmNC(benchmark)

  addpath(genpath(pwd),genpath('../utils/'));

  if nargin < 1
    benchmark = 'editable';
  end

  params = feval(benchmark);
  
  c4n = params.c4n;
  n4e = params.n4e;
  n4sDb = params.n4sDb;
  n4sNb = params.n4sNb;
  f = params.f;
  exactSolutionKnown = params.exactSolutionKnown; 
  % contains INITIAL data for the experiment
  
  % TODO make this work
  if params.showPlots
    figVisible = 'off';
  else
    figVisible = 'on';
  end
  % set(0,'DefaultFigureVisible','off');
  
  h = 2^(-params.initalRefinementLevel);
  
  %% MAIN AFEM LOOP

  eta4lvl = [];
  nrDoF4lvl = [];
  error4lvl = [];

  while( true )
    % SOLVE
    
    n4s = computeN4s(n4e);

    u = interpolationNC(f,c4n,n4e,n4s);
    
    du = computeGradientNC(c4n,n4e,u);
    varLambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
    varLambda(isinf(varLambda)) = 0;
    varLambda(isnan(varLambda)) = 0;
    
    %  varLambda = zeros(size(n4e,1),2);
    %  u = zeros(size(n4s,1),1); 
    tic;
    [u,corrVec,energyVec,nrDoF] = ...
      tvRegPrimalDual(params,c4n, n4e, n4sDb, n4sNb, u, varLambda,...
      params.epsStop, h, params.initalRefinementLevel);
    time = toc; 
    nrDoF4lvl(end+1) = nrDoF;
   
    % ESTIMATE
    eta4e = estimateError4e(u,f,c4n,n4e,n4sDb,n4sNb,1,1)
    eta4lvl(end+1) = sum(eta4e);
    disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDoF),...
      '; estimator = ',num2str(eta4lvl(end))]); 

    if
    error4lvl(end+1) = sqrt(sum(error4eCRL2(c4n,n4e,uExact,u)));

    saveResults('CR',expName,dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
      terminate,time,corrVec,energyVec,tau,miscMsg,nrDoF,...
      true,nrDoF4lvl,eta4lvl,error4lvl);

    % Check Termination
    if nrDoF >= minNrDoF, break, end;

    % MARK
    n4sMarked = markBulk(n4e,eta4e,0.5);

    % REFINE
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);

    if polygonGeometry
      temp=unique(n4sDb);
      c4n(temp,:)=c4n(temp,:)./repmat(sqrt(c4n(temp,1).^2+c4n(temp,2).^2),1,2);
    end
  end

  output = struct;

  output.eta4lvl = eta4lvl;
  output.nrDoF4lvl = nrDoF4lvl;
  output.error4lvl = error4lvl;
end
