function params = editable %#ok<*MSNU,FNDEF>
%% DOC
% Editable prototype for benchmark files.
% Execute program/startAlgorithmNC(benchmark).m (benchmark name of 
% benchmark file) to run algorithm.
% Execute program/computeExactEnergyBV(...) to approximate the exact energy.
%
% editable.m
% output: params - 'struct' with fields created dependent on the choice of the
%                  following parameters (documentation of parameters in the
%                  DOCUMENTATION OF PARAMETERS section below)

%% PARAMETERS
  % misc. parameters (will affect performance)
  showPlots              = false; 
  plotModeGrayscale      = false; 
    % not effective if showPlots == false
  showProgress           = true; 
  degree4Integrate       = 10; 
  plotGivenFunctions     = false;
  refinementLevel4Plots  = 5; % 11 is very close to the limit
    % not effective if plotGivenFunctions == false
  debugIfError           = false;

  % AFEM parameters
  geometry               = 'BigSquare'; %#ok<NASGU>                     
    % set automatically to 'Square' if useImage == true
  initialRefinementLevel = 4;
  parTheta               = 0.5;
  minNrDof               = 5e3;
  useProlongation        = true;
  parGamma               = 1;
  d                      = 2;
    % this should remain 2

  % algorithm parameters
  u0Mode         = 'zeros'; 
  initialEpsStop = 1e-4; 
  parTau         = 1; 
  maxIter        = 1e10;

  % experiment parameters
  useImage               = false;
  imageName              = 'f2bawgnSnr20cameraman.tif'; %#ok<NASGU> 
    % not effective if useImage == false
    % whiteSquare.tif, cameraman.tif
  parAlpha               = 1e0; 
  parBeta                = 1;
  rhsGradientKnown       = true;
    % set automatically to false if useImage == true
  exactSolutionKnown     = true; %#ok<NASGU>
    % set automatically to false if useImage == true
  useExactEnergy         = true; %#ok<NASGU>
    % set automatically to false if exactSolutionKnown == false
  exactEnergy            = -2.058034062391; %#ok<NASGU> % 6 significant digits
    % set automatically to NaN if exactSolutionKnown == false
    % not effective if useExactEnergy == false
                              
  % information about experiment for saving and documentation.
  expName                = 'OPTIMIZE';
  %dirInfoName            = sprintf('epsStop=%e', initialEpsStop);
  dirInfoName            = sprintf('%s', ...
    datestr(now, 'yy_mm_dd_HH_MM_SS'));

  % function handles (not effective if useImage == true)
  function y = rightHandSide(x)
    y =  f01(x, [parAlpha, parBeta]);
  end

  function y = gradientRightHandSide(x)
    % not effective if useExactEnergy == false
    y =  f01Gradient(x, [parAlpha, parBeta]);
  end

  function y = exactSolution(x)
    % not effective if exactSolutionKnown == false
    y = f01ExactSolution(x, parBeta);
  end
   
%% DOCUMENTATION OF PARAMETERS
% misc. parameters (will affect performance)
%   showPlots             - 'logical' with value 1 if plots must be shown
%                           during the iteration and 0 else
%   plotModeGrayscale     - 'logical' with value 1 if plots during the
%                            iteration must be grayscale plots with view from
%                            above onto the x-y plane and 0 else 
%   showProgress          - 'logical' with value 1 if information about the
%                           iteration must be printed during the iteration and
%                           0 else
%   degree4Integrate      - 'uint64' containing the algebraic degree of
%                            exactness for integrate from the AFEM package that
%                            must be used for calculations
%   plotGivenFunctions    - 'logical' with value 1 if a plot of the right-hand
%                           side f and, if given, of the exact solution must
%                           be saved and, if showPlots == true, shown
%   refinementLevel4Plots - 'uint64' containing the refinement level of the 
%                           mesh the right-hand side f and, if given, exact
%                           solution must be drawn on
%   debugIfError          - 'logical' with value 1 if MATLAB must enter debug
%                            mode if an error occurs during runtime and 0 else
% 
% AFEM parameters
%   geometry               - 'string'/'char array with exactly one row' 
%                            containing the name of the geometry the
%                            computation must be one
%   initialRefinementLevel - 'uint64' containing the initial refinement level
%                            the geometry must be loaded with by loadGeometry
%                            for level 0 of AFEM
%   parTheta               - 'double' containing the bulk parameter for marking
%                            (1 for uniform)
%   minNrDof               - 'uint64' containing the minimal number of degrees
%                            of freedom the last level of AFEM must have before
%                            terminating
%   useProlongation        - 'logical' with value 1 if the prolongation of the
%                             solution of a level of AFEM to the refined mesh 
%                             for the next level must be used as initial value
%                             for the iteration on the next level
%   parGamma               - 'double' containing the parameter \gamma for
%                            the refinement indicator
%   d                      - 'uint64' containing the dimension
%
% algorithm parameters
%   u0Mode         - 'char array with exactly one row' containing the choice
%                    for the inital iterate for iteration on level 0 and, if
%                    useProlongation == false, for the iterations on all levels
%     Options:
%                  'zeros': CR function with all coefficients equal to 0
%       'interpolationRhs': CR interpolation of the right-hand side f to the
%                           mesh on the level
%   initialEpsStop - 'double' containing the initial epsStop for relevant for 
%                    the iteration on the first level of AFEM
%   parTau         - 'double' containing the parameter \tau from the algorithm
%   maxIter        - 'uint64' TODO
%
% experiment parameters
%   useImage           - 'logical' with value 1 if an image given by imageName
%                        in the folder '../utils/functions/images/' must be
%                        used as right-hand side f for the experiment
%   imageName          - 'char array with exactly one row' containing the name
%                        of the image 
%   parAlpha           - 'double' containing the parameter \alpha from the
%                        problem
%   parbeta            - 'double' containing the parameter \beta from the
%                        problem
%   rhsGradientKnown   - 'logical' with value 1 if the gradient of the 
%                        given right-hand side is given and 0 else
%   exactSolutionKnown - 'logical' with value 1 if the exact solution to the
%                        continuous problem for the given right-hand side is
%                        given and 0 else
%   useExactEnergy     - 'logical' with value 1 if the exact energy given by
%                        exactEnergy must be used during runtime
%   exactEnergy        - 'double' containing the exact energy of the continuous
%                        problem for exact solution in H^1_0
%
% information about experiment for saving and documentation.
%   expName     - 'char array with exactly one row' containing the name for the
%                 folder in '../../results/nonconforming/' where the results of
%                 the experiment must be saved in
%   dirInfoName - 'char array with exactly one row' containing the name for the 
%                 folder in '../../results/nonconforming/expName/' where the 
%                 results of the experiment must be saved in

%% BUILD STRUCT
% advanced, not of interest for mere usage of the program
  if useImage
    geometry = 'Square'; %#ok<UNRCH>
    rhsGradientKnown = false; 
    exactSolutionKnown = false; 
    imageName =  sprintf('../utils/functions/images/%s', imageName);
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
  params.d = d;
  params.parGamma = parGamma;
  params.initialEpsStop = initialEpsStop;
  params.useProlongation = useProlongation;      
  params.rhsGradientKnown = rhsGradientKnown; 
  params.exactSolutionKnown = exactSolutionKnown; 

  if exactSolutionKnown 
    params.uExact = @(x) exactSolution(x); %#ok<UNRCH>
    params.useExactEnergy = useExactEnergy;
    params.exactEnergy = exactEnergy;
  else
    params.uExact = @(x) NaN; %#ok<UNRCH>
    params.useExactEnergy = false;
    params.exactEnergy = NaN;
  end

  params.parTau = parTau; 
  params.maxIter = maxIter; 

  params.parAlpha = parAlpha;
  params.parBeta = parBeta;

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
    params.f = image2function(imageName, parAlpha); %#ok<UNRCH>
  else
    params.f = @(x) rightHandSide(x); 
    if rhsGradientKnown 
      params.gradF = @(x) gradientRightHandSide(x);
    end
  end
end
