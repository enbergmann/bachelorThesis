function startAlgorithmCR(benchmark)
%% DOC
% Loads a benchmark and starts the corresponding experiment with the
% nonconforming algorithm.
%
% startAlgorithmCR.m
% input: benchmark - 'string'/'char array with exactly one row' containing the 
%                    name of the benchmark the user wants to use (optional
%                    parameter, default value is 'editable').

%% INIT
  % initialize paths and load benchmark
  cd(fileparts(which('startAlgorithmCR')));
    % change to the directory 'nonconforming' where startAlgorithmCR.m should
    % have been called from for correctness of relative filepaths during
    % runtime
  addpath(genpath(pwd), ...
    genpath('../utils/'), ...
    genpath('../conforming/plot/'));

  if nargin<1, benchmark = 'editable'; end

  params = feval(benchmark);
  params.benchmark = benchmark;
  if params.debugIfError, dbstop if error; end

  % extract necessary parameters from params
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
  useProlongation = params.useProlongation;
  rhsGradientKnown = params.rhsGradientKnown;
  useExactEnergy = params.useExactEnergy;
  exactEnergy = params.exactEnergy;
  u0Mode = params.u0Mode;
  
  % initialize outputLvl structs (structs for AFEM output)

  outputLvlHidden.hMax = [];
  outputLvlHidden.hMin = [];
  outputLvlHidden.sumL1NormOfJumps = [];
  
  lvl = 0;

  outputLvlInfo.lvl = lvl;
  nrDof = []; 
  outputLvlInfo.nrDof = nrDof; 
  outputLvlInfo.nrIterations = [];
  time = [];
  outputLvlInfo.time = time;

  outputLvlError.lvl = lvl;
  if exactSolutionKnown, outputLvlError.error4lvl = []; end 
  outputLvlError.eta = [];
  outputLvlError.etaVol = [];
  outputLvlError.etaJumps = [];

  outputLvlEnergy.lvl = lvl;
  outputLvlEnergy.energy = [];
  if useExactEnergy, outputLvlEnergy.diffDiscExactE = []; end
  if rhsGradientKnown 
    outputLvlEnergy.gleb = []; 
    if useExactEnergy, outputLvlEnergy.diffGlebExactE = []; end
    outputLvlEnergy.diffGlebDiscreteE = [];
  end
  
  % initialize currData (struct with parameters and data for the level)
  currData.c4n = c4n;
  currData.n4e = n4e;
  currData.n4sDb = n4sDb;
  currData.n4sNb = n4sNb;

  n4s = computeN4s(n4e);
  currData.n4s = n4s;

  nrSides = size(n4s, 1);
  currData.nrSides = nrSides;

  length4s = computeLength4s(c4n, n4s);
  currData.length4s = length4s;

  switch u0Mode
    case 'zeros', u0 = zeros(nrSides, 1);
    case 'interpolationRhs', u0 = interpolationCR(params, currData, f); 
  end

  % TODO here sth must be done when there are different possibilities for
  % epsStop
  currData.epsStop = params.initialEpsStop;

%% MAIN
  while(true)
    % initialize remaining current data
    hMax = max(length4s);
    hMin = min(length4s);
    currData.hMax = hMax;
    currData.hMin = hMin;
    outputLvlHidden.hMax(end+1, 1) = hMax;
    outputLvlHidden.hMin(end+1, 1) = hMin;

    currData.area4e = computeArea4e(c4n, n4e);

    currData.s4n = computeS4n(n4e);
    currData.s4e = computeS4e(n4e);

    currData.mid4e = computeMid4e(c4n, n4e); 
      % NOTE mid4e only for debugging with plotInitialTriangulation from utils
    currData.s4n = computeS4n(n4e, n4s);
    currData.e4s = computeE4s(n4e);

    currData.nrElems = size(n4e, 1);

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

    nrDof(end+1, 1) = length(dof);  %#ok<AGROW>
    currData.nrDof = nrDof(end);
    outputLvlInfo.nrDof = nrDof;

    [currData.int1RHS4e, currData.int2RHS4e, currData.int3RHS4e, ...
      currData.intRHS4s] = ...
      integralsWithF4e(params, currData);
      % needed here and in error estimate function

    % TODO
    % compute epsStop dependend on information given in benchmark
    % e.g. scaled with meshsize
    %
    % RIGHT NOW its just the initial epsStop as termination crit.!!!!

    % SOLVE (and save output information about the iteration)
    tic;
    [u, output.corrVec, energyVec, output.otherCorr] = ...
      solvePrimalDualFormulation(params, currData, u0, varLambda);
    outputLvlInfo.time(end+1, 1) = toc;
    outputLvlEnergy.energy(end+1, 1) = energyVec(end);
    output.energyVec = energyVec;
    output.u = u;
    outputLvlInfo.nrIterations(end+1, 1) = length(energyVec);
    if useExactEnergy 
      outputLvlEnergy.diffDiscExactE(end+1, 1) = ...
        abs(exactEnergy - energyVec(end)); 
    end

    output.normDiffRhsSolCrSquared4e = ...
      computeNormDiffRhsSolCrSquared4e(params, currData, output);

    % compute guaranteed lower energy bound
    if rhsGradientKnown
      glebCurr = computeGleb(params, currData, output);
      outputLvlEnergy.gleb(end+1, 1) = glebCurr;
      outputLvlEnergy.diffGlebDiscreteE(end+1, 1) = energyVec(end) - glebCurr;
      if useExactEnergy
        outputLvlEnergy.diffGlebExactE(end+1, 1) = exactEnergy - glebCurr;
      end
    end

    % ESTIMATE
    
    l1NormOfJump4s = computeL1NormOfJump4s(currData, output);
    currData.l1NormOfJump4s = l1NormOfJump4s;
    outputLvlHidden.sumL1NormOfJumps(end+1, 1) = sum(l1NormOfJump4s);

    [eta4e, etaVol4e, etaJumps4e] = ...
      estimateErrorCR4e(params, currData, output);
    
    % compute error
    % TODO implement flag for different errors
    if exactSolutionKnown
      outputLvlError.error4lvl(end+1, 1) = ...
        sqrt(sum(error4eCRL2(c4n, n4e, uExact, u)));
    end

    outputLvlError.eta(end+1, 1) = sum(eta4e);
    outputLvlError.etaVol(end+1, 1) = sum(etaVol4e);
    outputLvlError.etaJumps(end+1, 1) = sum(etaJumps4e);

    % display information about the level
    clc;
    disp(struct2table(outputLvlInfo));
    disp(struct2table(outputLvlError));
    disp(struct2table(outputLvlEnergy));

    % save results
    saveResultsCR(params, currData, ...
      outputLvlInfo, outputLvlError, outputLvlEnergy, outputLvlHidden, output);

    % check termination and update level
    if nrDof(end) >= minNrDof, break; end

    lvl(end+1, 1) = lvl(end) + 1; %#ok<AGROW>
    outputLvlInfo.lvl = lvl;
    outputLvlError.lvl = lvl;
    outputLvlEnergy.lvl = lvl;

    % MARK
    if parTheta==1, n4sMarked = markUniform(n4e); 
    else, n4sMarked = markBulk(n4e, eta4e, parTheta); end

    % REFINE (and prolongate solution to new mesh if needed)
    c4nOld = c4n;
    n4eOld = n4e;

    [c4n, n4e, n4sDb, n4sNb] = refineRGB(c4n, n4e, n4sDb, n4sNb, n4sMarked);

    if useProlongation
      % compute inital value for the iteration on the next level if
      % useProlongation (needs to be done before projection of the nodes on the
      % edges if polygonMesh, else getParentSide will not work)
      u0 = computeRefinementExtensionCR(c4nOld, n4eOld, c4n, n4e, u);
    end

    if polygonMesh
      % project boundary nodes of the refined mesh outward onto the unit circle
      % if polygonMesh
      temp = unique(n4sDb);
      c4n(temp, :) = ...
        c4n(temp, :)./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
    end
    
    % update first information in currData
    currData.c4n = c4n;
    currData.n4e = n4e;
    currData.n4sDb = n4sDb;
    currData.n4sNb = n4sNb;

    n4s = computeN4s(n4e);
    currData.n4s = n4s;

    nrSides = size(n4s, 1);
    currData.nrSides = nrSides;
  
    length4s = computeLength4s(c4n, n4s);
    currData.length4s = length4s;

    if ~useProlongation 
      % compute inital value for the iteration on the next level if not
      % useProlongation
      switch u0Mode
        case 'zeros', u0 = zeros(nrSides, 1);
        case 'interpolationRhs', u0 = interpolationCR(params, currData, f); 
      end
    end

    %TODO epsStop should probably be updated right here and now
  end
end
