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

%% INITIALIZATION
  addpath(genpath(pwd), genpath('../utils/'));

  if nargin < 1
    benchmark = 'editable';
  end

  params = feval(benchmark);
  
  % TODO this might have to be shortenend later, because some stuff just gets
  % copied to currGeom
  
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
  
  % initialize remaining parameters and struct with information dependend solely
  % on the current geometry
  eta4lvl = [];
  nrDof4lvl = [];
  error4lvl = [];


  currData = struct;

  currData.c4n = c4n;
  currData.n4e = n4e;
  currData.n4sDb = n4sDb;
  currData.n4sNb = n4sNb;


  n4s = computeN4s(n4e);
  currData.n4s = n4s;

  length4s = computeLength4s(c4n, n4s);
  currData.length4s = length4s;

  % TODO this could have a flag for different options
  u0 = interpolationCR(currData, f);

  currData.epsStop = params.epsStop;

%%%%%%%%%%%%%%%

%% MAIN AFEM LOOP
  while(true)
    
    currData.hMax = max(length4s);
    currData.area4e = computeArea4e(c4n,n4e);

    currData.s4n = computeS4n(n4e);

    s4e = computeS4e(n4e);
    currData.s4e = s4e;

    currData.nrElems = size(n4e, 1);
    currData.nrSides = max(max(s4e));

    currData.gradsCR4e = computeGradsCR4e(currData);
    gradCRu0 = gradientCR(currData, u0);

    % TODO could have an option for different inital lambda
    varLambda = gradCRu0./repmat(sqrt(sum(gradCRu0.^2,2)),1,2);
    varLambda(isinf(varLambda)) = 0;

    dof = computeDofCR(currData);
    currData.dof = dof;

    nrDof = length(dof);
    currData.nrDof = nrDof;
    nrDof4lvl(end+1) = nrDof;

    [currData.int1RHS4e, currData.int2RHS4e, currData.int3RHS4e] = ...
      computeIntegrals(currData, f, 200);
    %needed here and in error estimate function

    % TODO
    % compute epsStop dependend on information given in benchmark
    % e.g. scaled with meshsize
    %
    % RIGHT NOW its just the inital epsStop as termination crit.!!!!

    % SOLVE
    tic;
%TODO continue here
    [u, corrVec, energyVec] = ...
      tvRegPrimalDual(params, currData, u0, varLambda);
      % TODO might change name later
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

    %TODO consider prolongation and stuff
    u0 = interpolationCR(f,c4n,n4e,n4s);
    gradCRu = gradientCR(c4n,n4e,u);
    % TODO
    %n4s = computeN4s(n4e);
    %the mesh should be update here c4n, n4e, n4sDb, n4sNb
    %so should u and du, pretty much everything that was known before the
    %'while' loop began
    %TODO stuff below has to be done here
    n4s = computeN4s(n4e);
    currData.n4s = n4s;
  
    length4s = computeLength4s(c4n, n4s);
    currData.length4s = length4s;
  end

  output = struct;

  output.eta4lvl = eta4lvl;
  output.nrDof4lvl = nrDof4lvl;
  output.error4lvl = error4lvl;
end


%TODO REMEMBER MATLAB ist copyOnWrite, so try not to change structs in functions
%to avoid the struct being copied
%might be necessary to use 'clear' to delete data later (e.g. after saving c4n in current)

% TODO for all my functions (dont want to touch afem stuff) compute all
% necessary stuff (in particular that is dependend on geometry) before and pass
% it to the functions i.e. my mentality will be efficiency >> memory usage

%TODO save all current geometry stuff in a 'current' struct
%comments then should read sth like
%current: struct, must contain the following fields
%            c4n - [copy allready existing comments]
%            n4e -
%            ... -
% TODO might want a 'never quit modus' to termination (just do the algorithm until
% manual termination by the user)
