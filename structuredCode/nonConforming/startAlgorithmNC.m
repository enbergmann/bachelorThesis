function [params, output] = startAlgorithmNC(benchmark)

  %% INITIALIZATION

  addpath(genpath(pwd), genpath('../utils/'));

  if nargin < 1
    benchmark = 'editable';
  end

  params = feval(benchmark);
  
  % extract parameters from params
  
  showPlots = params.showPlots;
  initalRefinementLevel = params.initalRefinementLevel;
  c4n = params.c4n;
  n4e = params.n4e;
  n4sDb = params.n4sDb;
  n4sNb = params.n4sNb;
  f = params.f;
  epsStop = params.epsStop;
  exactSolutionKnown = params.exactSolutionKnown; 
  uExact = params.uExact;
  polygonMesh = params.polygonMesh;
  expName = params.expName;
  
  % initialize remaining parameters 
  
  eta4lvl = [];
  nrDof4lvl = [];
  error4lvl = [];

  %% MAIN AFEM LOOP
  
  while(true)
    % SOLVE
    
    n4s = computeN4s(n4e);

    %TODO
    u = interpolationNC(f,c4n,n4e,n4s);
    
    du = computeGradientNC(c4n,n4e,u);
    varLambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
    varLambda(isinf(varLambda)) = 0;
    varLambda(isnan(varLambda)) = 0;

    
    %  varLambda = zeros(size(n4e,1),2);
    %  u = zeros(size(n4s,1),1); 
    
    tic;
    [u,corrVec,energyVec,nrDof] = ...
      tvRegPrimalDual(params,c4n, n4e, n4sDb, n4sNb, u, varLambda,...
      epsStop, initalRefinementLevel);
    time = toc; 
    nrDof4lvl(end+1) = nrDof;
   
    % ESTIMATE
    eta4e = estimateError4e(u,f,c4n,n4e,n4sDb,n4sNb,1,1)
    eta4lvl(end+1) = sum(eta4e);
    disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDof),...
      '; estimator = ',num2str(eta4lvl(end))]); 

    if exactSolutionKnown
      error4lvl(end+1) = sqrt(sum(error4eCRL2(c4n,n4e,uExact,u)));
    end

    saveResults('CR',expName,dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
      terminate,time,corrVec,energyVec,tau,miscMsg,nrDof,...
      true,nrDof4lvl,eta4lvl,error4lvl);

    % Check Termination
    if nrDof >= minNrDof, break, end;

    % MARK
    n4sMarked = markBulk(n4e,eta4e,0.5);

    % REFINE
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,n4sNb,n4sMarked);

    if polygonMesh
      temp=unique(n4sDb);
      c4n(temp,:)=c4n(temp,:)./repmat(sqrt(c4n(temp,1).^2+c4n(temp,2).^2),1,2);
    end
  end

  output = struct;

  output.eta4lvl = eta4lvl;
  output.nrDof4lvl = nrDof4lvl;
  output.error4lvl = error4lvl;
end
