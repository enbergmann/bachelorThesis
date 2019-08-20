function [params, output] = startAlgorithmCR(benchmark)
% Loads a benchmark and starts the corresponding experiment with the
% nonconforming algorithm.
%
% startAlgorithmCR.m
% input:  benchmark - 'string'/'char array' containing the name of the
%                     benchmark the user wants to use (optional parameter,
%                     default value is 'editable').
%
% output: params    - 'struct' containing the parameters obtained from the
%                     benchmark.
%         output    - 'struct' containing the results of the experiment.

% TODO for all my functions (dont want to touch afem stuff) compute all
% necessary stuff (in particular that is dependend on geometry) before and pass
% it to the functions i.e. my mentality will be efficiency >> memory usage

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

  n4s = computeN4s(n4e);

  %TODO compute gradients4e before, use them to compute du e.g. and also
  %pass them to tvReg etc. Already did this, see computeGradientsNCnew vs 
  %gradientCR (basically, replace gradientCr by this new thing)
  u = interpolationCR(f,c4n,n4e,n4s);
  du = gradientCR(c4n,n4e,u);

  %% MAIN AFEM LOOP
  
  while(true)
    % SOLVE
    
    % varLambda = bsxfun(@rdivide,du,sqrt(sum(du.^2,2))); 
    % varLambda(isnan(varLambda)) = 0;
    varLambda = du./repmat(sqrt(sum(du.^2,2)),1,2);
    varLambda(isinf(varLambda)) = 0;
    %TODO compute necessary information for algorithm that has further use,
    %e.g. information about the mesh (nrDof), i.e. everything that is not 
    %just used in tvReg but also after that
    %
    %should probably be saved in a 'current' struct to minimize number
    %of input for tvReg
    
    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));
    dof = computeDof(n4e,nrSides,n4sDb,n4sNb);

    nrDof = length(dof);
    nrDof4lvl(end+1) = nrDof;

    % TODO
    % compute epsStop dependend on information given in benchmark
    % e.g. scaled with meshsize

    %TODO
    tic;
    % [u,corrVec,energyVec,nrDof] = ...
    %   tvRegPrimalDual(params,c4n, n4e, n4sDb, n4sNb, u, varLambda,...
    %   epsStop, initalRefinementLevel);
    output = ...
      tvRegPrimalDual(params, c4n, n4e, n4sDb, n4sNb, u, varLambda,...
      epsStop, initalRefinementLevel);
    time = toc; 
   
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
    
    n4s = computeN4s(n4e);

    u = interpolationCR(f,c4n,n4e,n4s);
    du = gradientCR(c4n,n4e,u);
    %TODO consider prolongation and stuff
    %n4s = computeN4s(n4e);
    %the mesh should be update here c4n, n4e, n4sDb, n4sNb
    %so should u and du, pretty much everything that was known before the
    %'while' loop began
  end

  output = struct;

  output.eta4lvl = eta4lvl;
  output.nrDof4lvl = nrDof4lvl;
  output.error4lvl = error4lvl;
end
