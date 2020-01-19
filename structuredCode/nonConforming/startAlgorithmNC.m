function [params, output] = startAlgorithmCR(benchmark)
% Loads a benchmark and starts the corresponding experiment with the
% nonconforming algorithm.
%
% startAlgorithmCR.m
% input:  benchmark - 'string'/'char array with exactly one row' containing the 
%                     name of the benchmark the user wants to use (optional
%                     parameter, default value is 'editable').
%
% TODO
% output: params    - 'struct' containing the parameters obtained from the
%                     benchmark.
%         output    - 'struct' containing the results of the experiment.

%% INITIALIZATION
  addpath(genpath(pwd), genpath('../utils/'));

  if nargin < 1
    benchmark = 'editable';
  end

  params = feval(benchmark);
  
  % extract parameters from params
  c4n = params.c4n;
  n4e = params.n4e;
  n4sDb = params.n4sDb;
  n4sNb = params.n4sNb;

  f = params.f;
  exactSolutionKnown = params.exactSolutionKnown; 
  uExact = params.uExact;
  polygonMesh = params.polygonMesh;
  minNrDof = params.minNrDof;
  parTheta = params.parTheta;
  
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

  % TODO this could have a flag for different options for initial u0
  % which would also have to be in benchmark (initialU0)
  u0 = interpolationCR(currData, f);

  currData.epsStop = params.epsStop;

%% MAIN AFEM LOOP
  while(true)
    
    currData.hMax = max(length4s);
    currData.area4e = computeArea4e(c4n, n4e);

    currData.s4n = computeS4n(n4e);

    s4e = computeS4e(n4e);
    currData.s4e = s4e;

    currData.nrElems = size(n4e, 1);
    currData.nrSides = max(max(s4e));

    currData.mid4e = computeMid4e(c4n, n4e);
    currData.s4n = computeS4n(n4e,n4s);
    currData.e4s = computeE4s(n4e);

    currData.gradsCR4e = computeGradsCR4e(currData);
    gradCRu0 = gradientCR(currData, u0);

    [currData.stiMaCR, currData.maMaCR] = computeFeMatricesCR(currData);

    % TODO could have an option for different initial lambda
    % which would also have to be in benchmark (initialVarLambda)
    varLambda = gradCRu0./repmat(sqrt(sum(gradCRu0.^2, 2)), 1, 2);
    varLambda(isinf(varLambda)) = 0; % this is prob. unnecessary, since only
                                     %0/0=NaN happens by definition
                                     %leave it just to be save? doesn't hurt
    varLambda(isnan(varLambda)) = 0;

    dof = computeDofCR(currData);
    currData.dof = dof;

    nrDof = length(dof);
    currData.nrDof = nrDof;
    nrDof4lvl(end+1) = nrDof;

    [currData.int1RHS4e, currData.int2RHS4e, currData.int3RHS4e, ...
      currData.intRHS4s] = ...
      integralsWithF4e(currData, f, 200);
    %needed here and in error estimate function

    % TODO
    % compute epsStop dependend on information given in benchmark
    % e.g. scaled with meshsize
    %
    % RIGHT NOW its just the initial epsStop as termination crit.!!!!

    % SOLVE
    tic;
    % TODO not done yet (subfunctions, documentation)
    [u, corrVec, energyVec] = ...
      solvePrimalDualFormulation(params, currData, u0, varLambda);
    time = toc; 

    % ESTIMATE

    %TODO still need to comment and some other stuff
    eta4e = estimateErrorCR4e(params, currData, u);
    eta4lvl(end+1) = sum(eta4e);

    disp(['nodes/dofs: ',num2str(size(c4n,1)),'/',num2str(nrDof),...
      '; estimator = ',num2str(eta4lvl(end))]); 

    % TODO implement flag for different errors
    if exactSolutionKnown
      error4lvl(end+1) = sqrt(sum(error4eCRL2(c4n,n4e,uExact,u)));
    end

    %%%%%%%%%%%%%%%%%%%%%TODO needs to be done
    saveResults('CR',expName,dirInfoName,figVisible,message,c4n,n4e,u,red,alpha,delta,...
      terminate,time,corrVec,energyVec,tau,miscMsg,nrDof,...
      true,nrDof4lvl,eta4lvl,error4lvl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check Termination
    if nrDof >= minNrDof, break, end;

    % MARK
    n4sMarked = markBulk(n4e, eta4e, parTheta);

    % REFINE
    [c4n, n4e, n4sDb, n4sNb] = refineRGB(c4n, n4e, n4sDb, n4sNb, n4sMarked);
    currData.c4n = c4n;
    currData.n4e = n4e;
    currData.n4sDb = n4sDb;
    currData.n4sNb = n4sNb;

    if polygonMesh
      temp = unique(n4sDb);
      c4n(temp, :) = ...
        c4n(temp, :)./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
    end
    
    n4s = computeN4s(n4e);
    currData.n4s = n4s;
  
    length4s = computeLength4s(c4n, n4s);
    currData.length4s = length4s;

    %TODO consider prolongation and stuff, choices for start values of u0
    u0 = interpolationCR(currData, f);

    %TODO epsStop should probably be updated right here and now
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

% TODO think about when do I want to write functions for stuff
%  - multiple use
%  - e.g. saveScreenshot stuff (to shorten code and not comment it since it's
%    not necesary at all for the code, so just put it somewhere where nobody
%    needs to see it)
