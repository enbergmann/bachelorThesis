function params = editable
  % Editable prototype for benchmark files.
  % Execute program/startAlgorithmNC.m to run algorithm.


  %%% PARAMETERS %%%

  % AFEM parameters
  geometry               = 'BigSquare';
  parTheta               = 0.5;  % bulk param. (1 for uniform)
  initalRefinementLevel  = 0;

  % algorithm parameters
  minNrDof               = 1e6;
  alpha4Estimate         = 1;
  beta4Estimate          = 1;   
  epsStop                = 1e-3;
  stopCrit               = ["Exact Error Difference", "weighted energy difference"];
  useProlongation        = false;
  exactSolutionKnown     = false;
  useExactEnergy         = false; % only effective if exactSolutionKnown == true
  %TODO compute exact energy here and use it e.g. in showProgress instead of
  %-2.0580......
  parTau                 = 1/2;

  % experiment parameters
  parAlpha               = 1;
  parBeta                = 1;
  errorNorm              = ["L2", "energy"]; % TODO list options (likewise for some other params (think))
                                              %strArr
  saveScreenshots        = 0; % save screenshots every saveScreenshots
                              % iterations during algorithm, e.g. for the case it doesn't finish
                              % 0 means no screenshots will be saved

  % misc. parameters (will affect performance)
  showPlots              = false; % Show plots during computation?
  showProgress           = true; % Print output during computation?

  % Information about experiment for saving and documentation.
  expName                = 'testForBenchmark';
  dirInfoName = datestr(now,'yy_mm_dd_HH_MM_SS');
  miscMsg                = sprintf(['this\nis\nan\nexample',...
                           '\non\nhow\nthis\ncould\nlook']);

  % necessary function handles
  function val = rightHandSide(x)
    val =  g(x,1,1);
  end

  function initalValue(x)
    val = 0;
  end

  function exactSolution(x)
    % can be ignored if exact solution is unknown
    val = gUexact(x,1,1);
  end


  %%% BUILD STRUCT %%%
  % should not be of interest for mere usage of the program
  
  params = struct;

  params.geometry = geometry;

  if strcmp(geometry, 'Polygon')
    [params.c4n, params.n4e, params.n4sDb, params.n4sNb] = ...
      computeGeometryPolygon(initalRefinementLevel);
    params.polygonMesh = true;
  else
    [params.c4n, params.n4e, params.n4sDb, params.n4sNb] = ...
      loadGeometry(geometry, initalRefinementLevel);
    params.polygonMesh = false;
  end

  params.parTheta = parTheta;
  params.initalRefinementLevel = initalRefinementLevel;

  params.minNrDof = minNrDof;
  params.alpha4Estimate = alpha4Estimate;
  params.beta4Estimate = beta4Estimate;
  params.epsStop = epsStop;
  params.stopCrit = stopCrit;              
  params.useProlongation = useProlongation;      
  params.exactSolutionKnown = exactSolutionKnown; 

  if exactSolutionKnown
    params.uExact = @(x) exactSolution(x);
  else
    params.uExact = @(x) NaN;
  end

  if exactSolutionKnown 
    params.useExactEnergy = useExactEnergy;
    %TODO probably compute it here
  end

  params.parTau = parTau; 

  params.parAlpha = parAlpha;
  params.parBeta = parBeta;
  params.errorNorm = errorNorm;             
  params.saveScreenshots = saveScreenshots;       

  params.showPlots = showPlots;              
  params.showProgress = showProgress;          

  params.expName = expName;
  params.dirInfoName = dirInfoName;
  params.miscMsg = miscMsg;

  params.f = @(x) rightHandSide(x);
  params.u0 = @(x) initalValue(x);
end
