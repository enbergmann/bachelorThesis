function params = editable %#ok<*MSNU,FNDEF>
%% DOC
% Editable prototype for benchmark files.
% Execute program/startAlgorithmNC(benchmark).m (benchmark name of 
% benchmark file) to run algorithm.
% Execute program/computeExactEnergyBV(...) to approximate the exact energy.

%% PARAMETERS
  % misc. parameters (will affect performance)
  showPlots              = false; 
  plotModeGrayscale      = false; 
    % not effective if showPlots == false
  showProgress           = true; 
  degree4Integrate       = 10; 
  plotGivenFunctions     = true;
  refinementLevel4Plots  = 4; % 11 is very close to the limit
    % not effective if plotGivenFunctions == false
  debugIfError           = false;

  % showPlots             - 'logical' with value 1 if plots must be shown
  %                         during the iteration and 0 else
  % plotModeGrayscale     - 'logical' with value 1 if plots during the
  %                          iteration must be grayscale plots with view from
  %                          above onto the x-y plane and 0 else 
  % showProgress          - 'logical' with value 1 if information about the
  %                         iteration must be printed during the iteration and
  %                         0 else
  % degree4Integrate      - 'uint64' containing the algebraic degree of
  %                          exactness for integrate from the AFEM package that
  %                          must be used for calculations
  % plotGivenFunctions    - 'logical' with value 1 if a plot of the right-hand
  %                         side f and, if given, of the exact solution must
  %                         be saved and, if showPlots == true, shown
  % refinementLevel4Plots - 'uint64' containing the refinement level of the 
  %                         mesh the right-hand side f and, if given, exact
  %                         solution must be drawn on
  % debugIfError          - 'logical' with value 1 if MATLAB must enter debug
  %                         mode if an error occurs during runtime and 0 else

  % AFEM parameters
  geometry               = 'BigSquare'; %#ok<NASGU>                     
    % set automatically to 'Square' if useImage == true
  initialRefinementLevel = 0;
  parTheta               = 0.5;
  minNrDof               = 5e3;
  useProlongation        = true;
  beta4Estimate          = 1;
  n4Estimate             = 2;
    % this should remain 2

  % geometry               - 'string'/'char array with exactly one row' 
  %                          containing the name of the geometry the
  %                          computation must be one
  % initialRefinementLevel - 'uint64' containing the initial refinement level
  %                          the geometry must be loaded with by loadGeometry
  %                          for level 0 of AFEM
  % parTheta               - 'double' containing the bulk parameter for marking
  %                          (1 for uniform)
  % minNrDof               - 'uint64' containing the minimal number of degrees
  %                          of freedom the last level of AFEM must have before
  %                          terminating
  % useProlongation        - 'logical' with value 1 if the prolongation of the
  %                           solution of a level of AFEM to the refined mesh 
  %                           for the next level must be used as initial value
  %                           for the iteration on the next level
  % beta4Estimate          - 'double' containing the parameter \beta from the
  %                          problem
  % n4Estimate             - 'uint64' containing the dimension

  % algorithm parameters
  u0Mode         = 'zeros'; 
  initialEpsStop = 1e-2; 
  stopCrit       = ["Exact Error Difference", ...
                    "weighted energy difference"]; 
  parTau         = 1/2;

  % u0Mode         - 'char array with exactly one row' containing the choice
  %                  for the inital iterate for iteration on level 0 and, if
  %                  useProlongation == false, for the iterations on all levels
  %   Options:
  %                'zeros': CR function with all coefficients equal to 0
  %     'interpolationRhs': CR interpolation of the right-hand side f to the
  %                         mesh on the level
  % initialEpsStop - 'double' containing the initial epsStop for relevant for 
  %                  the iteration on the first level of AFEM
  % stopCrit       - TODO
  % parTau         - 'double' containing the parameter \tau from the algorithm

  % experiment parameters
  useImage               = false;
  imageName              = ...
    '../utils/functions/images/whiteSquare.tif'; %#ok<NASGU> 
    % whiteSquare, cameraman
    % TODO append path automatically
  addNoise               = false; %#ok<NASGU>
    % TODO probably add ability to denoise rhs and images (might need case
    % distinction by considering useImage flag)
    % TODO noise type (see MATLAB imnoise)
  blurWidth              = 1; %#ok<NASGU>
    % TODO comment
    % 1 for nothing  TODO think about it, it might change stuff
  parAlpha               = 1e0; %1e4 for image example 
   % TODO why does the analytic example is broken for 1e4
  parBeta                = 1;
  exactSolutionKnown     = true; %#ok<NASGU>
  useExactEnergy         = true; %#ok<NASGU>
    % only effective if exactSolutionKnown == true

		% just write it in from the file per hand, with like 10 digits or 
		% sth.. Think about it.
    % TODO for now this is also the flag to use gradient of f to compute a
    % guaranteed lower energy bound
  % TODO how should the exactEnergy be written into here
  exactEnergy            = -2.05805109; %#ok<NASGU>
    % four significant digits
                              % TODO 
  saveScreenshots        = 0; % save screenshots every saveScreenshots
                              % iterations during algorithm, e.g. for the case
                              % it doesn't finish
                              % 0 means no screenshots will be saved

  % Information about experiment for saving and documentation.
  expName                = 'commentCommon';
  dirInfoName            = sprintf('%s', ...
    datestr(now, 'yy_mm_dd_HH_MM_SS'));
  errorNorm              = ["L2", "energy"]; 
    %TODO list options (likewise for some other params (think)) strArr maybe
    %altErrorNorm and always compute L2 error because it's useful for the
    %estimator (or always L2 and H1 error, if useful, but they probably only
    %differ bei some constant hence one of those is sufficient (check this))
    %TODO might call this altError and give only ine option (i.e. there's
    %always L2 error and maybe altError) UNTIL THEN error4lvl is a fine name
    %and is always L2 
    %TODO maybe allow only a fixed amounts of different errors, like only two,
    %so one can choose L2 and H1 for example so sth like error4lvl,
    %errorAlt4lvl
    %TODO probably always have the error for which the estimator is an upper
    %bound

  % function handles (can be ignored if useImage)
  function y = rightHandSide(x)
    % TODO pasted-graphic-2.tiff does have a calculation formula to calculate f
    % from some given function u(r) --> other examples possible (easier even?)
    % y =  middleSquare(x);
    y =  f01(x, [parAlpha, parBeta]);
  end

  function y = gradientRightHandSide(x)
    y =  f01Gradient(x, [parAlpha, parBeta]);
  end

  function y = exactSolution(x)
    % can be ignored if exactSolutionKnown == false
    y = f01ExactSolution(x, parAlpha);
  end

%% BUILD STRUCT
% should not be of interest for mere usage of the program
  if useImage
    geometry = 'Square'; %#ok<UNRCH>
      % TODO for now only on square possible, but 
      % might be possible on rectangles in general (useful 
      % though?)
    exactSolutionKnown = false; % NOTE there won't be exact solutions for
                                % real images
  end

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
  params.n4Estimate = n4Estimate;
  params.beta4Estimate = beta4Estimate;
  params.initialEpsStop = initialEpsStop;
  params.stopCrit = stopCrit;              
  params.useProlongation = useProlongation;      
  params.exactSolutionKnown = exactSolutionKnown; 

  if exactSolutionKnown 
    params.uExact = @(x) exactSolution(x); %#ok<UNRCH>
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
    params.uExact = @(x) NaN; %#ok<UNRCH>
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
  if showPlots, params.figVisible = 'on'; %#ok<UNRCH>
  else, params.figVisible = 'off'; end
  params.plotModeGrayscale = plotModeGrayscale;
  params.showProgress = showProgress;          
  params.plotGivenFunctions = plotGivenFunctions;
  params.refinementLevel4Plots = refinementLevel4Plots;
  params.debugIfError = debugIfError;

  params.expName = expName;
  params.dirInfoName = dirInfoName;

  params.u0Mode = u0Mode;
  if useImage
    params.f = image2function(imageName, parAlpha, ...
      addNoise, blurWidth); %#ok<UNRCH>
    %TODO rewrite and change name
  else
    noise = 0; %#ok<NASGU,UNRCH> 
    %if addNoise, noise = .4; end;
    params.f = @(x) rightHandSide(x); %+ noise*(rand-.5);
    %params.f = @(x) rightHandSide(x) + noise*(rand(size(x,1))-.5);
    %   but this is a different noise on every level
    %  TODO what if add noise, this obv doesnt work since it only adds one 
    %  random number
    if useExactEnergy %TODO this is also possible without exact energy
                      % just needs to rewrite some stuff in saveResults etc.
      params.gradF = @(x) gradientRightHandSide(x);
    end
  end
end
