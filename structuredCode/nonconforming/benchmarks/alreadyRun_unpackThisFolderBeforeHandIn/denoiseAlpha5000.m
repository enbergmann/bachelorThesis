function params = editable %#ok<*MSNU,FNDEF>
%% DOC
% Editable prototype for benchmark files.
% Execute startAlgorithmCR(benchmark) (benchmark name of 
% benchmark file) from folder ./nonconforming/ to run algorithm.
% Use ./nonconforming/computeExactEnergyBV.m to approximate the exact
% energy of a given problem with known exact solution.
%
% editable.m
% output: params - 'struct' with fields created dependent on the choice of the
%                  following parameters (documentation of parameters in the
%                  DOCUMENTATION OF PARAMETERS section below)

%% PARAMETERS
  % misc. parameters (will affect performance)
  showPlots              = false; 
  plotModeGrayscale      = false; % not effective if showPlots==false
  showProgress           = true; 
  degree4Integrate       = 10; 
  plotGivenFunctions     = false;
  refinementLevel4Plots  = 8; % 11 is very close to the limit
    % not effective if plotGivenFunctions==false
  debugIfError           = false;

  % AFEM parameters
  geometry               = 'BigSquare'; %#ok<NASGU>                     
    % set automatically to 'Square' if useImage==true
  initialRefinementLevel = 0;
  parTheta               = 0.5;
  minNrDof               = 1e8;
  useProlongation        = true;
  parGamma               = 1;
  d                      = 2; % this should remain 2

  % algorithm parameters
  u0Mode  = 'zeros'; %'interpolationInSi'; 'zeros';
  epsStop = 1e-4; 
  parTau  = 1; 
  maxIter = 1e12;

  % experiment parameters
  useImage               = true;
  imageName              = 'f2bawgnSnr15cameraman.tif'; %#ok<NASGU> 
    % not effective if useImage==false
  parAlpha               = 5e3; 
  parBeta                = 1;
  inSiGradientKnown      = false;
    % set automatically to false if useImage==true
  exactSolutionKnown     = false; %#ok<NASGU>
    % set automatically to false if useImage==true
  useExactEnergy         = false; %#ok<NASGU>
    % set automatically to false if exactSolutionKnown==false
  exactEnergy            = -2.058034062391; %#ok<NASGU> % 6 significant digits
    % set automatically to NaN if exactSolutionKnown==false
    % not effective if useExactEnergy==false
                              
  % information about experiment for saving and documentation
  expName                = 'denoiseSNR15';
  dirInfoName            = sprintf('parAlpha=%.1e', parAlpha);

  % function handles (not effective if useImage==true)
  function y = inputSignal(x)
    y =  f01(x, [parAlpha, parBeta]);
  end

  function y = gradientInputSignal(x)
    % not effective if inSiGradientKnown==false
    y =  f01Gradient(x, [parAlpha, parBeta]);
  end

  function y = exactSolution(x)
    % not effective if exactSolutionKnown==false
    y = f01ExactSolution(x, parBeta);
  end
   
%% DOCUMENTATION OF PARAMETERS
% misc. parameters (will affect performance)
%   showPlots             - 'logical' with value 1 if plots must be shown
%                           during the iteration and 0 else
%   plotModeGrayscale     - 'logical' with value 1 if plots during the
%                           iteration must be grayscale plots with view from
%                           above onto the xy-plane and 0 else 
%   showProgress          - 'logical' with value 1 if information about the
%                           iteration must be printed during the iteration and
%                           0 else
%   degree4Integrate      - 'uint64' containing the algebraic degree of
%                           exactness for integrate from the AFEM package that
%                           must be used for calculations
%   plotGivenFunctions    - 'logical' with value 1 if a plot of the input
%                           signal f and, if given, of the exact solution must
%                           be saved and, if showPlots==true, shown
%   refinementLevel4Plots - 'uint64' containing the refinement level of the 
%                           mesh the input signal f and, if given, exact
%                           solution must be drawn on
%   debugIfError          - 'logical' with value 1 if MATLAB must enter debug
%                           mode if an error occurs during runtime and 0 else
% 
% AFEM parameters
%   geometry               - 'string'/'char array with exactly one row' 
%                            containing the name of the geometry the
%                            computation must be one
%   initialRefinementLevel - 'uint64' containing the initial refinement level
%                            of the geometry loaded by loadGeometry for level 0
%                            of AFEM
%   parTheta               - 'double' containing the bulk parameter for marking
%                            (1 for uniform)
%   minNrDof               - 'uint64' containing the minimal number of degrees
%                            of freedom the last level of AFEM must have before
%                            terminating
%   useProlongation        - 'logical' with value 1 if a prolongation of the
%                            solution of a level of AFEM to the refined mesh
%                            for the next level must be used as initial value
%                            for the iteration on the next level
%   parGamma               - 'double' containing the parameter \gamma for
%                            the refinement indicator
%   d                      - 'uint64' containing the dimension
%
% algorithm parameters
%   u0Mode  - 'char array with exactly one row' containing the choice for the
%             inital iterate for the iteration on level 0 and, if
%             useProlongation==false, for the iterations on all levels
%     Options:
%                   'zeros': CR function with all coefficients equal to 0
%       'interpolationInSi': CR interpolation of the input signal f to the mesh
%                            on the level
%   epsStop - 'double' containing \epsilon_{stop} for the termination criterion
%             of the iteration
%   parTau  - 'double' containing the parameter \tau for the iteration
%   maxIter - 'uint64' containing the number of iteration steps the iteration
%             must not exceed
%
% experiment parameters
%   useImage           - 'logical' with value 1 if an image given by imageName
%                        in the folder '../utils/functions/images/' must be
%                        used as input signal f for the experiment
%   imageName          - 'char array with exactly one row' containing the name
%                        of the image 
%   parAlpha           - 'double' containing the parameter \alpha from the
%                        problem
%   parbeta            - 'double' containing the parameter \beta from the
%                        problem (if existent)
%   inSiGradientKnown  - 'logical' with value 1 if the gradient of the given
%                        input signal is given and 0 else
%   exactSolutionKnown - 'logical' with value 1 if the exact solution to the
%                        continuous problem for the given input signal is given
%                        and 0 else
%   useExactEnergy     - 'logical' with value 1 if the exact energy given by
%                        exactEnergy must be used during runtime
%   exactEnergy        - 'double' containing the exact energy of the continuous
%                        problem for exact solution in H^1_0
%
% information about experiment for saving and documentation.
%   expName     - 'char array with exactly one row' containing the name for the
%                 folder in '../../results/nonconforming/' where the results of
%                 runs of the experiment must be saved in
%   dirInfoName - 'char array with exactly one row' containing the name for the 
%                 folder in '../../results/nonconforming/expName/' where the
%                 results of a particular run of the experiment must be saved
%                 in

%% BUILD STRUCT
% advanced, not of interest for mere usage of the program
  if useImage
    geometry = 'Square'; %#ok<UNRCH>
    inSiGradientKnown = false; 
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
  params.epsStop = epsStop;
  params.useProlongation = useProlongation;      
  params.inSiGradientKnown = inSiGradientKnown; 
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
    params.f = @(x) inputSignal(x); 
    if inSiGradientKnown 
      params.gradF = @(x) gradientInputSignal(x);
    end
  end
end
