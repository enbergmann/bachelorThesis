function  [u, corrVec, energyVec, otherCorr] = ...
    solvePrimalDualFormulation(params, currData, u, varLambda) 
%% DOC
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
%                        saveScreenshots: 'uint64' containing the information
%                                         that every saveScreenshots-th
%                                         iteration a screenshot of the results
%                                         of the iteration must be saved and 0
%                                         if no screenshots must be saved
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
%         corrVec   - '(1 x nrIterations) - dimensional double array' where the
%                     j-th column contains the corr for the j-th iteration step
%         energyVec - '(1 x nrIterations) - dimensional double array' where the
%                     j-th column contains the discrete energy of the j-th
%                     iterate 

%% INIT
  % extract necessary parameters from params
  parTau = params.parTau;
  parAlpha = params.parAlpha;
  showProgress = params.showProgress;
  exactEnergy = params.exactEnergy;
  useExactEnergy = params.useExactEnergy;
  saveScreenshots = params.saveScreenshots;
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

  A = stiMaCR/parTau+parAlpha*maMaCR; 

  gradCRu = gradientCR(currData, u);
  E = computeDiscreteEnergyCR(params, currData, u, gradCRu);

  v = zeros(nrSides, 1);    

  corrVec = [];
  energyVec = [E];

  % prepare comparison between different termination criteria
  otherCorr.eNcAbsDiffVec = [];

  hMin = currData.hMin;
  C = maMaCR + hMin*stiMaCR;
  otherCorr.bar15TerminationVec = [];
  otherCorr.bar15TerminationWithoutL2Vec = []; 
    % p. 314, Alg. 10.1; p. 316, Prop. 10.8; p. 317, Section 10.2.4
  otherCorr.bar12TerminationVec = []; 
    % p. 1173, Section 6.2
  otherCorr.bar12TerminationSqrtVec = []; 
    
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
  while true
    % compute basic information for the next itertion step
    gradCRv = gradientCR(currData, v);
    M = varLambda + parTau*(gradCRu + parTau*gradCRv);

    varLambda = M./repmat(max(1, sqrt(sum(M.^2, 2))), 1, 2);
    varLambda(isinf(varLambda)) = 0; % this is prob. unnecessary, since only
                                     %0/0=NaN happens by definition
                                     %leave it just to be save? doesn't hurt
    varLambda(isnan(varLambda)) = 0;
    
    % compute right-hand side

    % for k = 1, ..., nrSides
    %   b(k) = ...
    %     (temp, \nabla_{NC}\psi_k)_{L^2(\Omega)} ...
    %     + intRHS4s(k)
    % with
    %   temp = \nabla_{NC}u_{j-1} - \Lambda_j    and
    %   intRHS4s(k) = (f, \psi_k)_L^2(\Omega)
    % furthermore, utilizing that temp and \nabla\psi_k are constant on any
    % triangle and therefore also temp\cdot \nabla\psi_k,
    %   (temp, \nabla_{NC}\psi_k)_{L^2(\Omega)} ...
    %     = \int_\Omega temp\cdot\nabla_{NC}\psi_k dx 
    %     = \sum{T\in\mathcal{T}} \int_T temp\cdot\nabla\psi_k dx
    %     = \sum{T\in\mathcal{T}} temp\cdot\nabla\psi_k \int_T 1 dx
    %     = \sum{T\in\mathcal{T}} temp\cdot\nabla\psi_k |T|
    % and since psi_k = 0 on \mathcal{T}\setminus{T_{+}, T_{-}}
    %   (temp, \nabla_{NC}\psi_k)_{L^2(\Omega)} ...
    %     = temp\cdot\nabla\psi_k|_{T_{+}} |T_{+}| ...
    %       + temp\cdot\nabla\psi_k|_{T_{-}} |T_{-}|
    % for an inner edge E_k = T_{+}\cap T_{-} and
    %   (temp, \nabla_{NC}\psi_k)_{L^2(\Omega)} ...
    %     = temp\cdot\nabla\psi_k|_{T_{+}} |T_{+}|
    
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
      % TODO here case distinction for termination criteria needs to be done
    corr = sqrt(dtU'*stiMaCR*dtU); % Only gradients
      % TODO gradCRu already computed, so using the stiMa might be unnecessary
      % (cf. Tiens code)
      
    % for comparision of termination criteria

    otherCorr.eNcAbsDiffVec(end+1) = abs(E-ENew);
    otherCorr.bar15TerminationVec(end+1) = sqrt(dtU'*C*dtU);
    otherCorr.bar15TerminationWithoutL2Vec(end+1) = ...
      sqrt(hMin*dtU'*stiMaCR*dtU); 
    otherCorr.bar12TerminationVec(end+1) = dtU'*maMaCR*dtU/energyVec(1); 
    otherCorr.bar12TerminationSqrtVec(end+1) = ...
      sqrt(dtU'*maMaCR*dtU)/energyVec(1); 

    % update data
    u = uNew;
    E = ENew;
    energyVec(end+1) = E; %#ok<AGROW>
    corrVec(end+1) = corr; %#ok<AGROW>

    % show miscellaneous information
    if showProgress
      fprintf(repmat('\b', 1, lineLength));
      lineLength = fprintf('%e      %f        %d', ...
        corr, E, length(corrVec));
    end

    if saveScreenshots > 0 && mod(length(energyVec), saveScreenshots) == 0
      % TODO this function is not written yet, do it next time it's needed
      saveScreenshot();
    end

    if showPlots
      clf('reset');
      if plotModeGrayscale, plotGrayscale(c4n, n4e, mean(u(s4e), 2));
      else, plotCR(c4n,n4e,uNew); end
    end

    % check termination
    if corr<epsStop, break; end
  end
end
