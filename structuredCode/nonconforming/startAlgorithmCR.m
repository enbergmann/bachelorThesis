function startAlgorithmCR(benchmark)
%% DOC
% Loads a benchmark and starts the corresponding experiment with the
% nonconforming algorithm.
%
% startAlgorithmCR.m
% input: benchmark - 'string'/'char array with exactly one row' containing the 
%                    name of the benchmark the user wants to use (optional
%                    parameter, default value is 'editable')

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
  inSiGradientKnown = params.inSiGradientKnown;
  useExactEnergy = params.useExactEnergy;
  exactEnergy = params.exactEnergy;
  u0Mode = params.u0Mode;
  
  % initialize outputLvl structs (structs for AFEM output)
  outputLvlHidden.hMin = [];
  outputLvlHidden.hMax = [];
  outputLvlHidden.sumL1NormOfJumps = [];
  if exactSolutionKnown, outputLvlHidden.normDiffExactSolJ1DiscSol = []; end
  outputLvlHidden.normDiffDiscSolJ1DiscSol = [];
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
  outputLvlEnergy.gueb = []; 
  if inSiGradientKnown 
    outputLvlEnergy.gleb = []; 
    outputLvlEnergy.diffGuebGleb = []; 
    if useExactEnergy, outputLvlEnergy.diffExacEGleb = []; end
    outputLvlHidden.diffDiscEGleb = [];
  end
  if useExactEnergy, outputLvlEnergy.absDiffDiscExacE = []; end
  
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
    case 'interpolationInSi', u0 = interpolationCR(params, currData, f); 
  end

%% MAIN
  while(true)
    % initialize remaining current data
    hMin = min(length4s);
    hMax = max(length4s);
    currData.hMin = hMin;
    currData.hMax = hMax;
    outputLvlHidden.hMin(end+1, 1) = hMin;
    outputLvlHidden.hMax(end+1, 1) = hMax;

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

    [currData.stiMaCR, maMaCR] = computeFeMatricesCR(currData);
    currData.maMaCR = maMaCR;

    varLambda = gradCRu0./repmat(sqrt(sum(gradCRu0.^2, 2)), 1, 2);
    varLambda(isinf(varLambda)) = 0; % this is prob. unnecessary, since only
                                     % 0/0 = NaN happens by definition
                                     % remains here just to be save
    varLambda(isnan(varLambda)) = 0;

    dof = computeDofCR(currData);
    currData.dof = dof;

    nrDof(end+1, 1) = length(dof);  %#ok<AGROW>
    currData.nrDof = nrDof(end);
    outputLvlInfo.nrDof = nrDof;

    [currData.int1InSi4e, currData.int2InSi4e, currData.int3InSi4e, ...
      currData.intInSi4s] = ...
      integralsWithInSi(params, currData);
      % needed here and in error estimate function

    % SOLVE (and save output information about the iteration)
    tic;
    [u, outputLvlInfo.nrIterations(end+1, 1), ...
      output.corrVec, energyVec, output.otherCorr] = ...
      solvePrimalDualFormulation(params, currData, u0, varLambda);
    outputLvlInfo.time(end+1, 1) = toc;
    outputLvlEnergy.energy(end+1, 1) = energyVec(end);
    output.energyVec = energyVec;
    output.u = u;
    if useExactEnergy 
      outputLvlEnergy.absDiffDiscExacE(end+1, 1) = ...
        abs(exactEnergy - energyVec(end)); 
    end

    output.normDiffInSiSolCrSquared4e = ...
      computeNormDiffInSiSolCrSquared4e(params, currData, output);

    % compute guaranteed energy bounds
    uJ1 = computeJ1(n4e, n4sDb, u);
    uJ1 = courant2CR(n4e, uJ1);
    uJ1GradCR = gradientCR(currData, uJ1);
    guebCurr = computeDiscreteEnergyCR(params, currData, uJ1, uJ1GradCR);
    outputLvlEnergy.gueb(end+1, 1) = guebCurr;
    if inSiGradientKnown
      glebCurr = computeGleb(params, currData, output);
      outputLvlEnergy.gleb(end+1, 1) = glebCurr;
      outputLvlEnergy.diffGuebGleb(end+1, 1) = guebCurr - glebCurr; 
      outputLvlHidden.diffDiscEGleb(end+1, 1) = energyVec(end) - glebCurr;
      if useExactEnergy
        outputLvlEnergy.diffExacEGleb(end+1, 1) = exactEnergy - glebCurr;
      end
    end
    temp = u - uJ1;
    outputLvlHidden.normDiffDiscSolJ1DiscSol(end+1, 1) = ...
      sqrt(temp'*maMaCR*temp);

    % ESTIMATE
    
    l1NormOfJump4s = computeL1NormOfJump4s(currData, output);
    currData.l1NormOfJump4s = l1NormOfJump4s;
    outputLvlHidden.sumL1NormOfJumps(end+1, 1) = sum(l1NormOfJump4s);

    [eta4e, etaVol4e, etaJumps4e] = ...
      estimateErrorCR4e(params, currData, output);
    
    % compute error
    if exactSolutionKnown
      outputLvlError.error4lvl(end+1, 1) = ...
        sqrt(sum(error4eCRL2(c4n, n4e, uExact, u)));
      outputLvlHidden.normDiffExactSolJ1DiscSol(end+1, 1) = ...
        sqrt(sum(error4eCRL2(c4n, n4e, uExact, uJ1)));
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
    n4sDbOld = n4sDb;

    [c4n, n4e, n4sDb, n4sNb] = refineRGB(c4n, n4e, n4sDb, n4sNb, n4sMarked);

    if useProlongation
      % compute inital value for the iteration on the next level if
      % useProlongation (needs to be done before projection of the nodes on the
      % edges if polygonMesh, else getParentSide will not work)
      u0 = prolongationJ1(c4nOld, n4eOld, n4sDbOld, c4n, n4e, u);
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
        case 'interpolationInSi', u0 = interpolationCR(params, currData, f); 
      end
    end
  end
end
