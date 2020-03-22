function  [u,corrVec,energyVec] = ...
    solvePrimalDualFormulation(params, currData, u, varLambda) 
  % TODO think about how the output should be handled, then think about 
  % documentation

% Computes the piecewise gradient of the Crouzeix-Raviart function v (vanishing
% in the midpoints of boundary edges) with respect to the triangulation given
% by [c4n, n4e].
%
% gradientCR.m
% input:  params   - 'struct' with fields:
%         currData - 'struct' with fields:
%                         nrElems: number of elements
%                             s4e: sides for elements
%                       gradsCR4e: gradients of side based Crouzeix-Raviart
%                                  basis functions for all elements
%         v        - 'function_handle' of the function whose piecewise gradient 
%                     is to be computed
%
% output: u  - '(number of sides x 2)-dimensional double array' where 
%                     the j-th row contains the gradient of v on the j-th
%                     triangle 

  % extract necessary parameters from params
  parTau = params.parTau;
  parAlpha = params.parAlpha;
  showProgress = params.showProgress;
  exactEnergy = params.exactEnergy;
  useExactEnergy = params.useExactEnergy;
  saveScreenshots = params.saveScreenshots;
  showPlots = params.showPlots;
  plotModeGrayscale = params.plotModeGrayscale;

  % extract necessary data from currData
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

  % INITIALIZATION 
  
  % extract necessary data

  % initialize remaing parameters

  % firstScreenshot = datestr(now, 'yy_mm_dd_HH_MM_SS');

  A = stiMaCR/parTau+parAlpha*maMaCR; %TODO here could be an h in front of stiMaCR
  %C = maMaCR + h*stiMaCR;

  gradCRu = gradientCR(currData, u);

  v = zeros(nrSides, 1);    

  corr = epsStop+1; 
  corrVec = [];
  energyVec = [];

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

  while corr > epsStop
    gradCRv = gradientCR(currData, v);
    M = varLambda + parTau*(gradCRu + parTau*gradCRv);

    varLambda = M./repmat(max(1, sqrt(sum(M.^2, 2))), 1, 2);
    varLambda(isinf(varLambda)) = 0;
    
    % compute RHS
    
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

    %% Solve System
    uNew = zeros(nrSides, 1);
    uNew(dof) = A(dof, dof)\b(dof);
    v = (uNew - u)/parTau;        

    %% Check Termination
    gradCRu = gradientCR(currData, uNew);
    ENew = computeDiscreteEnergyCR(params, currData, uNew, gradCRu);

    % TODO here case distinction for termination criteria needs to be done
    dtU = (u - uNew)/parTau; 
    %corr = sqrt(dtU'*C*dtU); % Bartels termination criterion
    corr = sqrt(dtU'*stiMaCR*dtU); % Only gradients
      % TODO gradCRu already computed, so using the stiMa might be unnecessary
      % (cf. Tiens code)

    u = uNew;
    E = ENew;
    energyVec(end+1) = E; %#ok<AGROW>
    corrVec(end+1) = corr; %#ok<AGROW>

    if showProgress
      % TODO use structs and disp table maybe or just table (there must be
      % a getTable and replace feature)

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
      % subplot(1, 2, 1)
      % cla('reset');
      % plotCR(c4n, n4e, uNew);
      % subplot(1, 2, 2)
      % cla('reset');
      % plotGrayscale(c4n, n4e, mean(u(s4e), 2));
    end
  end
end
