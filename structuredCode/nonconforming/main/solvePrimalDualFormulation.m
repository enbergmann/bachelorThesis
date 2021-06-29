function  [u, nrIter, corrVec, energyVec, otherCorr] = ...
    solvePrimalDualFormulation(params, currData, u, varLambda) 
%% DOC TODO update
% Executes the nonconfomring iteration on the triangulation given
% by [c4n, n4e] with initial values u and varLambda.
%
% solvePrimalDualFormulation.m
% input: params    - 'struct' with fields:
%                                 parTau: 'double' containing the parameter 
%                                         \tau from the problem
%                               parAlpha: 'double' containing the parameter 
%                                         \alpha from the problem
%                           showProgress: 'logical' with value 1 if progress
%                                         must be printed during the iteration
%                                         and 0 else
%                            exactEnergy: 'double' containing the exact energy
%                                         of the minimizer of the continuous
%                                         problem
%                         useExactEnergy: 'logical' with value 0 if exactEnergy
%                                         must be ignored and 1 else
%                              showPlots: 'logical' with value 1 if plots
%                                         must be shown during the iteration
%                                         and 0 else
%                      plotModeGrayscale: 'logical' with value 1 if shown plots
%                                         must be plotted as grayscale plots
%                                         with view from above onto the x-y
%                                         plane
%        currData  - 'struct' with fields:
%                        stiMaCR: '(nrSides x nrSides)-dimensional sparse
%                                 double array' where the k-th entry in the
%                                 j-th row is the L2-scalar product of the
%                                 piecewise gradient of the CR-basis function
%                                 w.r.t. the k-th edge with the piecewise
%                                 gradient of the CR-basis function w.r.t. the
%                                 j-th edge
%                         maMaCR: '(nrSides x nrSides)-dimensional sparse
%                                 double array' where the k-th entry in the
%                                 j-th row is the L2-scalar product of the
%                                 CR-basis function w.r.t. the k-th edge with
%                                 the CR-basis function w.r.t. the j-th edge
%                            c4n: coordinates for nodes
%                            n4e: nodes for elements  
%                        nrSides: number of sides
%                        epsStop: 'double' containing the minimal value corr
%                                 must reach for the iteration to terminate
%                      gradsCR4e: '(3 x 2 x nrElems)-dimensional double array'
%                                 where the j-th row of the k-th matrix
%                                 contains the gradient of the local CR-basis
%                                 function w.r.t. the j-th edge of the k-th
%                                 element
%                            s4e: sides for elements
%                         area4e: areas for elements
%                       intRHS4s: '(nrSides x 1)-dimensional double array'
%                                 where the j-th entry is the integral of f
%                                 times the CR-basis function w.r.t. the j-th
%                                 edge
%                            dof: '(1 x nrDof)-dimensional double array' where
%                                 the j-th column contains the number of the
%                                 j-th degree of freedom
%                          nrDof: number of degrees of freedom
%                            e4s: elements for sides
%        u         - '(nrSides x 1)-dimensional double array' where the j-th
%                    row contains the CR coefficient of the initial value u for
%                    the iteration
%        varLambda - '(nrELems x 2)-dimensional double array' where the j-th
%                    row contains the initial value for Lambda for the
%                    iteration on the j-th element of the triangulation
%
% output: u         - '(nrSides x 1)-dimensional double array' where the j-th
%                     row contains the CR coefficient of the final iterate of
%                     the iteration, i.e. the approximated solution
%         corrVec   - '(1 x nrIter) - dimensional double array' where the
%                     j-th column contains the corr for the j-th iteration step
%         energyVec - '(1 x nrIter) - dimensional double array' where the
%                     j-th column contains the discrete energy of the j-th
%                     iterate 

%% INIT
  % extract necessary parameters from params
  parTau = params.parTau;
  maxIter = params.maxIter;
  parAlpha = params.parAlpha;
  showProgress = params.showProgress;
  exactEnergy = params.exactEnergy;
  useExactEnergy = params.useExactEnergy;
  showPlots = params.showPlots;
  plotModeGrayscale = params.plotModeGrayscale;


  % extract necessary information from currData
  stiMaCR = currData.stiMaCR;
  maMaCR = currData.maMaCR;
  c4n = currData.c4n;
  n4e = currData.n4e;
  nrSides = currData.nrSides;
  epsStop = currData.epsStop;
  gradsCR4e = currData.gradsCR4e;  
  s4e = currData.s4e;
  area4e = currData.area4e;
  intRHS4s = currData.intRHS4s;
  dof = currData.dof;
  nrDof = currData.nrDof;
  e4s = currData.e4s;

  % initialize further variables
  % firstScreenshot = datestr(now, 'yy_mm_dd_HH_MM_SS');

  chunkSize = 10000; % size of preallocated entries for dynamic arrays
  A = stiMaCR/parTau+parAlpha*maMaCR; 

  gradCRu = gradientCR(currData, u);
  E = computeDiscreteEnergyCR(params, currData, u, gradCRu);

  v = zeros(nrSides, 1);    

  corrVec = NaN([1, chunkSize]);
  energyVec = [E, NaN([1, chunkSize])];

  % prepare comparison between different termination criteria
  otherCorr.eNcAbsDiffVec = NaN([1, chunkSize]);

  hMin = currData.hMin;
  C = maMaCR + hMin*stiMaCR;
  otherCorr.bar15TerminationVec = NaN([1, chunkSize]);
  otherCorr.bar15TerminationWithoutL2Vec = NaN([1, chunkSize]); 
    % p. 314, Alg. 10.1; p. 316, Prop. 10.8; p. 317, Section 10.2.4
  otherCorr.bar12TerminationSqrtVec = NaN([1, chunkSize]); 
    % p. 1173, Section 6.2
  
  % prepare computation of b
  elemPlusForSides = e4s(:, 1);
    % the j-th component contains the number of T_{+} for the j-th edge of the
    % mesh
  elemMinusForInnerSides = e4s(dof, 2);
    % the j-th component contains the number of T_{-} for the j-th inner edge
    % of the mesh
  
  sidesPlus = s4e(elemPlusForSides, :);
    % the k-th entry of the j-th row contains the global number of the k-th
    % local side of T_{+} for the j-th edge of the mesh
  sidesMinus = s4e(elemMinusForInnerSides, :);
    % the k-th entry of the j-th row contains the global number of the k-th
    % local side of T_{-} for the j-th inner edge of the mesh

  gradCrOnElemPlusForSides = zeros(nrSides, 2);
  for side = 1:nrSides
    localNr = sidesPlus(side, :) == side;
      % local number of side on T_{+}
    gradCrOnElemPlusForSides(side, :) = ...
      gradsCR4e(localNr, :, elemPlusForSides(side));
      % the j-th row contains the gradient of the j-th CR basis function w.r.t.
      % the j-th edge of the mesh on T_{+}
  end

  gradCrOnElemMinusForInnerSides = zeros(nrDof, 2); 
  for nrInnerSide = 1:nrDof 
    localNr = sidesMinus(nrInnerSide, :) == dof(nrInnerSide);
      % local number of the nrInnerSide-th side on T_{-}
    gradCrOnElemMinusForInnerSides(nrInnerSide, :) = ...
      gradsCR4e(localNr, :, elemMinusForInnerSides(nrInnerSide));
      % the j-th row contains the gradient of the CR basis function w.r.t. the
      % j-th inner edge of the mesh on T_{-}
  end

  % start printing progress and initialize figure if showPlots
  if showProgress
    fprintf('\n========================================\n\n');
    fprintf('Current iteration on a mesh with \n\n');
    fprintf('nrDof: %d\n', nrDof);
    fprintf('epsStop: %e\n', epsStop);
    if useExactEnergy, fprintf('Exact energy: %f\n', exactEnergy); end
    fprintf('\n========================================\n\n');
    fprintf('    corr           energy         step\n');
    fprintf('   ------         --------       ------ \n');
    lineLength = 0;
  end
  if showPlots, figure; end

%% MAIN
  nrIter = 0;
  while true
    nrIter = nrIter + 1;

    % compute basic information for the next itertion step
    gradCRv = gradientCR(currData, v);
    M = varLambda + parTau*(gradCRu + parTau*gradCRv);

    varLambda = M./repmat(max(1, sqrt(sum(M.^2, 2))), 1, 2);
    varLambda(isinf(varLambda)) = 0; % this is prob. unnecessary, since only
                                     %0/0=NaN happens by definition
                                     %leave it just to be save? doesn't hurt
    varLambda(isnan(varLambda)) = 0;
    
    % compute right-hand side
    temp4e = (gradCRu/parTau - varLambda);
    b = intRHS4s + ...
      area4e(elemPlusForSides).*sum(...
      temp4e(elemPlusForSides, :).*gradCrOnElemPlusForSides, 2);
    b(dof) = b(dof) + ...
      area4e(elemMinusForInnerSides).*sum(...
      temp4e(elemMinusForInnerSides, :).*gradCrOnElemMinusForInnerSides, 2);

    % solve system
    uNew = zeros(nrSides, 1);
    uNew(dof) = A(dof, dof)\b(dof);
    v = (uNew - u)/parTau;        

    % compute information dependent on the new iterate uNew
    gradCRu = gradientCR(currData, uNew);
    ENew = computeDiscreteEnergyCR(params, currData, uNew, gradCRu);

    dtU = (u - uNew)/parTau; % = -v;
    corr = sqrt(dtU'*stiMaCR*dtU); % Only gradients
      
    % check if new memory for dynamic arrays needs to be preallocated
    if nrIter > length(corrVec)
      corrVec = [corrVec, NaN([1, chunkSize])]; %#ok<AGROW>
      energyVec = [energyVec, NaN([1, chunkSize])]; %#ok<AGROW>
      otherCorr.eNcAbsDiffVec = ...
        [otherCorr.eNcAbsDiffVec, NaN([1, chunkSize])]; 
      otherCorr.bar15TerminationVec = ...
        [otherCorr.bar15TerminationVec, NaN([1, chunkSize])]; 
      otherCorr.bar15TerminationWithoutL2Vec = ...
        [otherCorr.bar15TerminationWithoutL2Vec, ...
        NaN([1, chunkSize])]; 
      otherCorr.bar12TerminationSqrtVec = ...
        [otherCorr.bar12TerminationSqrtVec, NaN([1, chunkSize])]; 
    end

    % compute other possible termination criteria for comparison
    otherCorr.eNcAbsDiffVec(nrIter) = abs(E-ENew); 
    otherCorr.bar15TerminationVec(nrIter) = sqrt(dtU'*C*dtU);
    otherCorr.bar15TerminationWithoutL2Vec(nrIter) = ...
      sqrt(hMin*dtU'*stiMaCR*dtU);
    otherCorr.bar12TerminationSqrtVec(nrIter) = ...
      sqrt(dtU'*maMaCR*dtU); 

    % update data
    u = uNew;
    E = ENew;
    energyVec(nrIter + 1) = E; 
    corrVec(nrIter) = corr;

    % show miscellaneous information
    if showProgress
      fprintf(repmat('\b', 1, lineLength));
      lineLength = fprintf('%e      %f        %d', ...
        corr, E, nrIter);
    end

    if showPlots
      clf('reset');
      if plotModeGrayscale, plotGrayscale(c4n, n4e, mean(u(s4e), 2));
      else, plotCR(c4n,n4e,uNew); end
    end

    % check termination
    if corr < epsStop || nrIter >= maxIter, break; end
  end
  % cut down dynamic arrays to their actual size
  corrVec = corrVec(1:nrIter);
  energyVec = energyVec(1:nrIter + 1);
  otherCorr.eNcAbsDiffVec = otherCorr.eNcAbsDiffVec(1:nrIter);
  otherCorr.bar15TerminationVec = ...
    otherCorr.bar15TerminationVec(1:nrIter);
  otherCorr.bar15TerminationWithoutL2Vec = ...
    otherCorr.bar15TerminationWithoutL2Vec(1:nrIter);
  otherCorr.bar12TerminationSqrtVec = ...
    otherCorr.bar12TerminationSqrtVec(1:nrIter);
end
