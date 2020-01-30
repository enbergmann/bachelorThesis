function params = editable
% Editable prototype for benchmark files.
% Execute program/startAlgorithmNC.m to run algorithm.
% Execute program/computeExactEnergyBV(...) to approximate the exact energy.

  % TODO useProlongation = true

%% PARAMETERS

  % AFEM parameters
  geometry               = 'Polygon'; % TODO does this work with adaptivity?
                                      % sth goes terribly wrong
  parTheta               = 0.5;  % bulk param. (1 for uniform)
  initialRefinementLevel = 0;

  % algorithm parameters
  minNrDof               = 1e3;
  alpha4Estimate         = 1;
  beta4Estimate          = 1;   
  epsStop                = 1e-2;
  stopCrit               = ["Exact Error Difference", ...
                            "weighted energy difference"];
  useProlongation        = true; % TODO
  exactSolutionKnown     = true;
  useExactEnergy         = true; % only effective if exactSolutionKnown == true
			   % just write it in from the file per hand, with like 10 digits or 
			   % sth.. Think about it.
  parTau                 = 1/2;

  % experiment parameters
  parAlpha               = 1;
  parBeta                = 1;
  errorNorm              = ["L2", "energy"]; % TODO list options (likewise for
                                             % some other params (think))
                                              %strArr
                              % TODO 
  saveScreenshots        = 0; % save screenshots every saveScreenshots
                              % iterations during algorithm, e.g. for the case it doesn't finish
                              % 0 means no screenshots will be saved

  % misc. parameters (will affect performance)
  degree4Integrate       = 20; % Algebraic degree of exactness for integrate
                               % from the AFEM package
  showPlots              = true; % Show plots during computation?
  showProgress           = true; % Print output during computation?

  % Information about experiment for saving and documentation.
  expName                = 'testForBenchmark';
  dirInfoName            = datestr(now, 'yy_mm_dd_HH_MM_SS');
  miscMsg                = sprintf(['this\nis\nan\nexample', ...
                                    '\non\nhow\nthis\ncould\nlook']);

  % TODO rename the called funtions at some time to more suitable names
  % necessary function handles
  function val = rightHandSide(x)
    val =  g(x, [1,1]);
  end

  function val = initalValue(x)
    val = 0;
  end

  function val = exactSolution(x)
    % can be ignored if exactSolutionKnown == false
    val = gUexact(x, [1,1]);
  end

  % TODO how should the exactEnergy be written into the programm
  exactEnergy = -2.05805109; % four significant digits

%% BUILD STRUCT
% should not be of interest for mere usage of the program
  params = struct;

  params.geometry = geometry;

  if strcmp(geometry, 'Polygon')
    [params.c4n, params.n4e, params.n4sDb, params.n4sNb] = ...
      computeGeometryPolygon(initialRefinementLevel);
    params.polygonMesh = true;
  else
    [params.c4n, params.n4e, params.n4sDb, params.n4sNb] = ...
      loadGeometry(geometry, initialRefinementLevel);
    params.polygonMesh = false;
  end

  params.parTheta = parTheta;
  params.initialRefinementLevel = initialRefinementLevel;

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
    params.exactEnergy = exactEnergy;
    % just compute it extrenalyy and only set the flag to true, if the energy is
    % given (read really exact valye from file)
    % TODO continue here, when the function is great


    %TODO probably compute it here
    %TODO compute exact energy here and use it e.g. in showProgress instead of
    % drop useExactEnergy and just use exactEnergy = 3 type params
    % use insan(exactEnergy) to find out whether exact energy is know or not, 
    % i.e. exactEnergy = NaN, if it is not known
    %
    % function that coputes exact energy should simply refine the mesh mesh
    % uniformly until it has a million dofs or something and then compute the 
    % energy and use it for here (might take some time, but only once, hence 
    % who cares)
    % this functions must prob. be written without struct s.t. it can be called
    % from here
    %
    % OR
    % compute exact energy in results where it is needed, but where
    % meshes are already computed to a certain degree (con: might still need
    % exact energy during programm (maybe in a termination criterion, even
    % though this might be stupid since there likely won't be results guaranteen 
    % sth))
    %-2.0580.....
    %
    %TODO
    %best solution: first time on mesh and given solution compute it very well
    %and save it (compute until some timer runs out, with message, saying that
    %it's computed at the moment)
    %
    %use struct [field = geometry, structInner] with structInner [field
    % = exactSolutionName,
    %intArray] with 1x2 intArray [nrDof, Error]
    %
    %use significant decimals at termination criterion
  else
    params.useExactEnergy = false;
    params.exactEnergy = NaN;
  end

  params.parTau = parTau; 

  params.parAlpha = parAlpha;
  params.parBeta = parBeta;
  params.errorNorm = errorNorm;             
  params.saveScreenshots = saveScreenshots;       

  params.degree4Integrate = degree4Integrate;
  params.showPlots = showPlots;              
  if showPlots
    params.figVisible = 'on';
  else
    params.figVisible = 'off';
  end
  params.showProgress = showProgress;          

  params.expName = expName;
  params.dirInfoName = dirInfoName;
  params.miscMsg = miscMsg;

  params.f = @(x) rightHandSide(x);
  params.u0 = @(x) initalValue(x);
end
