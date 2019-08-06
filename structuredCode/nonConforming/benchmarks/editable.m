function params = editable
  % Editable prototype for benchmark files.
  % Execute program/startAlgorithmNC.m to run algorithm.


  %%% PARAMETERS %%%

  % AFEM parameters
  geometry               = 'BigSquare';
  parTheta                  = 0.5;  % bulk param. (1 for uniform)
  initalRefinementLevel  = 0;

  % algorithm parameters
  minNrDoF               = 1e6;
  epsStop                = 1e-3;
  stopCrit               = 'Exact Error Difference';
  useProlongation        = false;
  exactSolutionKnown     = false;
  parTau                 = 1/2;

  % experiment parameters
  parAlph                = 1;
  parBeta                = 1;
  errorNorm              = 'L2'; % TODO list options (likewise for some other params (think))
  saveScreenshots        = 0; % save screenshots every saveScreenshots
                              % iterations during algorithm, e.g. for the case it doesn't finish
                              % 0 means no screenshots will be saved

  % misc. parameters (can affect performance)
  showPlots              = false; % Show plots during computation?
  showProgress           = true; % Print output during computation?

  % Information about experiment for saving and documentation.
  expName = 'testForBenchmark'
  dirInfoName = datestr(now,'yy_mm_dd_HH_MM_SS');
  miscMsg = sprintf('this\nis\nan\nexample');
  miscMsg = sprintf('%s\non\nhow\nthis\ncould\nlook');

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
  params.theta = theta
  params.initalRefinementLevel = initalRefinementLevel  
  params.f = @(x) rightHandSide(x);
  params.u0 = @(x) initalValue(x);
  if exactSolutionKnown
    params.uExact = @(x) exactSolution(x);
  else
    params.uExact = @(x) NaN;
  end

  if strcmp(geometry, 'Polygon')
    [params.c4n, params.n4e, params.n4sDb, params.n4sNb] = ...
      computeGeometryPolygon(initalRefinementLevel);
    params.polygonMesh = true;
  else
    [params.c4n, params.n4e, params.n4sDb, params.n4sNb] = ...
      loadGeometry(geometry, initalRefinementLevel);
    params.polygonMesh = false;
  end

  params.miscMsg = miscMsg;
end
