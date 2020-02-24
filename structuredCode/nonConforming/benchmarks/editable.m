function params = editable
% Editable prototype for benchmark files.
% Execute program/startAlgorithmNC.m to run algorithm.
% Execute program/computeExactEnergyBV(...) to approximate the exact energy.

% TODO document (or make clear) dependencies, i.e. which flags are
%      automatically set by which other flags
% TODO computeGleb flag, not effective for image probably due to missing 
%      gradients. For H^1_0 examples there always will be gradients probably
%      on the other hand.

%% PARAMETERS

  % misc. parameters (will affect performance)
  showPlots              = false; 
    % Show plots during iteration?
  showProgress           = true; 
    % Print output during iteration?
  plotModeGrayscale      = false; 
    % Use plotGrayscale instead of plotCR during iteration? Only effective if
    % showProgress.
  degree4Integrate       = 20; 
    % algebraic degree of exactness for integrate from the AFEM package

  % AFEM parameters
  geometry               = 'BigSquare'; 
    % not necessary if useImage (for now)
  parTheta               = 0.5;  
    % bulk param. (1 for uniform)
  initialRefinementLevel = 0;
  minNrDof               = 1e4;
  useProlongation        = true; 
  beta4Estimate          = 1;   

  % algorithm parameters
  epsStop                = 1e-2; % TODO sth about updating it depending on
                                 %      mesh size
  stopCrit               = ["Exact Error Difference", ...
                            "weighted energy difference"];
  parTau                 = 1/2;

  % experiment parameters
  useImage               = false;
  imageName              = '../utils/functions/images/cameraman.tif'; 
  addNoise               = false;
    % TODO probably add ability to denoise rhs and images (might need case
    % distinction by considering useImage flag)
    % TODO noise type (see MATLAB imnoise)
  parAlpha               = 1e0; %1e4 for image example 
   % TODO why does the analytic example is broken for 1e4
  parBeta                = 1;
  exactSolutionKnown     = true;
  useExactEnergy         = true; % only effective if exactSolutionKnown == true
		% just write it in from the file per hand, with like 10 digits or 
		% sth.. Think about it.
    % TODO for now this is also the flag to use gradient of f to compute a
    % guaranteed lower energy bound
  % TODO how should the exactEnergy be written into here
  exactEnergy            = -2.05805109; % four significant digits
  errorNorm              = ["L2", "energy"]; % TODO list options (likewise for
                                             % some other params (think))
                                              %strArr
                                            %maybe altErrorNorm and always
                                            %compute L2 error because it's 
                                            %useful for the estimator (or
                                            %always L2 and H1 error, if 
                                            %useful, but they probably only 
                                            %differ bei some constant hence
                                            %one of those is sufficient (check
                                            %this))
                              % TODO 
  saveScreenshots        = 0; % save screenshots every saveScreenshots
                              % iterations during algorithm, e.g. for the case it doesn't finish
                              % 0 means no screenshots will be saved

  % Information about experiment for saving and documentation.
  expName                = 'testDenoising';
  dirInfoName            = datestr(now, 'yy_mm_dd_HH_MM_SS');
  miscMsg                = sprintf(['this\nis\nan\nexample', ...
                                    '\non\nhow\nthis\ncould\nlook']);


  % function handles 
  function y = initalValue(x)
    y = 0;
  end

  % function handles that can be ignored if useImage
  % TODO pasted-graphic-2.tiff does have a calculation formula to calculate
  % f from some given function u(r) --> other examples possible (easier even?)
  function y = rightHandSide(x)
    y =  f01(x, [parAlpha, parBeta]);
  end

  function y = gradientRightHandSide(x)
    y =  f01Gradient(x, [parAlpha, parBeta]);
  end

  function y = exactSolution(x)
    % can be ignored if exactSolutionKnown == false
    y = f01ExactSolution(x, [parAlpha]);
  end


%% BUILD STRUCT
% should not be of interest for mere usage of the program
  if useImage
    geometry = 'Square'; % TODO for now only on square possible, but 
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
  params.beta4Estimate = beta4Estimate;
  params.epsStop = epsStop;
  params.stopCrit = stopCrit;              
  params.useProlongation = useProlongation;      
  params.exactSolutionKnown = exactSolutionKnown; 

  if exactSolutionKnown 
    params.uExact = @(x) exactSolution(x);
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
    params.uExact = @(x) NaN;
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
  if showPlots, params.figVisible = 'on';
  else, params.figVisible = 'off'; end
  params.plotModeGrayscale = plotModeGrayscale;
  params.showProgress = showProgress;          

  params.expName = expName;
  params.dirInfoName = dirInfoName;
  params.miscMsg = miscMsg;

  params.u0 = @(x) initalValue(x);
  if useImage
    params.f = rhsImg(imageName, parAlpha, addNoise); 
    %TODO rewrite and change name
  else
    noise = 0;
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
