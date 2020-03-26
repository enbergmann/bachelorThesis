function  [u, corrVec, energyVec] = ...
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
%                        nrElems: number of elements
%                            s4e: sides for elements
%                         area4e: areas for elements
%                       intRHS4s: '(nrSides x 1)-dimensional double array'
%                                 where the j-th entry is the integral of f
%                                 times the CR-basis function w.r.t. the j-th
%                                 edge
%                            dof: '(1 x nrDof)-dimensional double array' where
%                                 the j-th column contains the number of the
%                                 j-th degree of freedom
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
  nrElems = currData.nrElems;
  s4e = currData.s4e;
  area4e = currData.area4e;
  intRHS4s = currData.intRHS4s;
  dof = currData.dof;

  % initialize further variables
  % firstScreenshot = datestr(now, 'yy_mm_dd_HH_MM_SS');

  A = stiMaCR/parTau+parAlpha*maMaCR; 
    %TODO here could be an h in front of stiMaCR
  %C = maMaCR + h*stiMaCR;

  gradCRu = gradientCR(currData, u);

  v = zeros(nrSides, 1);    

  corrVec = [];
  energyVec = [];

  % start printing progress and initialize figure if showPlots
  if showProgress
    fprintf('\n========================================\n\n');
    fprintf('Current iteration on a mesh with \n\n');
    fprintf('nrDof: %d\n', length(dof));
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
    % b = zeros(nrSides, 1);

    % % NOTE 1st iteration of the code
    % for elem = 1 : nrElems
    %   bLocal = (gradCRu(elem, :)/parTau - varLambda(elem, :))...
    %     *gradsCR4e(:, :, elem)';  
    %   b(s4e(elem, :)) = b(s4e(elem, :)) + area4e(elem)*bLocal'; % right-hand side
    % end

    % NOTE 2nd iteration of the code
    %bLocal = zeros(nrElems, 3);
    % for elem = 1 : nrElems
    %   bLocal(elem, :) = bTemp(elem, :)*gradsCR4e(:, :, elem)';  
    % end

    % TODO rewrite it better, 
    % thought: sum works for depth dimension af array, so one could add first
    % and then multiply (or not, dunno)
    bTemp = (gradCRu/parTau - varLambda); % bTemp4e
    bRe = reshape(bTemp', 2*nrElems, 1);
    gRe = reshape(permute(gradsCR4e, [2 3 1]), 2*nrElems, 3);
    bLocal = area4e.*reshape(sum(reshape(bRe.*gRe, 2, nrElems*3)), nrElems, 3);

    % % NOTE 2nd iteration stuff
    % for elem = 1 : nrElems
    %   b(s4e(elem, :)) = b(s4e(elem, :)) + bLocal(elem, :)'; 
    %     % right-hand side
    % end
    % b = b + intRHS4s;

    % one side can only exist 1 or 2 times, not more (use unique fist and last)
    % this solution is highly CR0 dependend, because it needs dof and inner
    % edges to be the exact same (more general one would need to compute the
    % inner edges and use them instead of dof)

    [s4eSorted, s4eSortedInd] = sort(s4e(:));
    [~, s4eSortedUniqueIndFirst] = unique(s4eSorted);
    [~, s4eSortedUniqueIndLast] = unique(s4eSorted, 'last');
    s4eSortedUniqueFirst = s4eSortedInd(s4eSortedUniqueIndFirst);
    s4eSortedUniqueLast = s4eSortedInd(s4eSortedUniqueIndLast);
    bLocalFirst = bLocal(s4eSortedUniqueFirst);
    bLocalLast = bLocal(s4eSortedUniqueLast);
    b = intRHS4s + bLocalFirst;
    b(dof) = b(dof) + bLocalLast(dof);
      % without dof outer edges (non-dof for CR0) would be counted two times

    % solve system
    uNew = zeros(nrSides, 1);
    uNew(dof) = A(dof, dof)\b(dof);
    v = (uNew - u)/parTau;        

    % compute information dependent on the new iterate uNew
    gradCRu = gradientCR(currData, uNew);
    ENew = computeDiscreteEnergyCR(params, currData, uNew, gradCRu);

    dtU = (u - uNew)/parTau; 
      % TODO here case distinction for termination criteria needs to be done
    %corr = sqrt(dtU'*C*dtU); % Bartels termination criterion
    corr = sqrt(dtU'*stiMaCR*dtU); % Only gradients
      % TODO gradCRu already computed, so using the stiMa might be unnecessary
      % (cf. Tiens code)

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
