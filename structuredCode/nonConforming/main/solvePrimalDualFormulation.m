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
  saveScreenshots = params.saveScreenshots;
  showPlots = params.showPlots;
  imageGiven = params.imageGiven;

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
  firstScreenshot = datestr(now, 'yy_mm_dd_HH_MM_SS');

  A = stiMaCR/parTau+parAlpha*maMaCR; %TODO here could be an h in front of stiMaCR
  %C = maMaCR + h*stiMaCR;

  gradCRu = gradientCR(currData, u);

  v = zeros(nrSides, 1);    

  corr = epsStop+1; 
  corrVec = [];
  energyVec = [];
  E = 1;

  while corr > epsStop
    gradCRv = gradientCR(currData, v);
    M = varLambda + parTau*(gradCRu + parTau*gradCRv);

    varLambda = M./repmat(max(1, sqrt(sum(M.^2, 2))), 1, 2);
    varLambda(isinf(varLambda)) = 0;
    
    % compute RHS
    
    % TODO recalculate this at some point, Tiens looks different (smarter prob.)
    b = zeros(nrSides, 1);
    for elem = 1 : nrElems
      bLocal = (gradCRu(elem, :)/parTau - ...
        varLambda(elem, :))*gradsCR4e(:, :, elem)';  
      b(s4e(elem, :)) = b(s4e(elem, :)) + area4e(elem)*bLocal'; % right-hand side
    end
    b = b + intRHS4s;

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

    if showProgress
      % TODO make this prettier (also number of iteration)
      % use structs and disp table
      fprintf('corr/epsStop: %e / %e\n', corr, epsStop);
      format long; % TODO change that in the fprintf %.8g or something
      fprintf('E = %f, E_exact = %f\n', E, exactEnergy);
      format short;
      fprintf('============================== \n');
    end

    u = uNew;
    E = ENew;
    energyVec(end+1) = E;
    corrVec(end+1) = corr;

    if saveScreenshots > 0 & mod(length(energyVec), saveScreenshots) == 0
      % TODO this function is not written yet, do it next time it's needed
      saveScreenshot();
    end

    if showPlots
      if imageGiven
        colormap gray;
        axis off;
        axis equal;
        % plotGreyScale(mean(u(s4e), 2), c4n, n4e)
        X1 = c4n(n4e(:, 1), 1);
        X2 = c4n(n4e(:, 2), 1);
        X3 = c4n(n4e(:, 3), 1);
        Y1 = c4n(n4e(:, 1), 2);
        Y2 = c4n(n4e(:, 2), 2);
        Y3 = c4n(n4e(:, 3), 2);
        X = [X1'; X2'; X3'];
        Y = [Y1'; Y2'; Y3'];
        patch(X,Y,mean(u(s4e), 2)','EdgeColor','none');
        drawnow
      else
        plotCR(c4n,n4e,uNew);
      end
      clf('reset');
      fprintf('\n')
    end
  end
end
