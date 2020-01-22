function computeExactEnergyBV(geometry, fStr, fStrParams, uStr, uStrParams, ...
    gradUStr, gradUStrParams, parAlpha, ...,
    minNrDof, minPrecision, degree4Integrate)
% Computes and saves an approximation (at least up to precision minPrecision
% and minNrDof dofs) of the exact BV energy of a function u, whose pointwise
% gradient is known, for a right-hand side f on a mesh given by geometry.
%
% computeExactEnergyBV.m
% input:  geometry         - 'char array with exactly one row' containing the 
%                            name of the geometry the user wants to approximate
%                            the exact energy on
%         fStr             - 'string'/'char array with exactly one row' 
%                            containing the name of the right-hand side f
%         fStrParams       - 'double array with exactly one row' containing the 
%                            necessary parameters to produce a function handle
%                            of f (0x0 is possible if f needs no further
%                            parameters)
%         uStr             - 'string'/'char array with exactly one row' 
%                            containing the name of the function u whose exact
%                            BV energy is to be approximated 
%         uStrParams       - 'double array with exactly one row' containing the 
%                            necessary parameters to produce a function handle
%                            of u (0x0 is possible if u needs no further
%                            parameters)
%         gradUStr         - 'string'/'char array with exactly one row' 
%                            containing the name of the gradient of the
%                            function u whose exact BV energy is to be
%                            approximated
%         gradUStrParams   - 'double array with exactly one row' containing the 
%                            necessary parameters to produce a function handle 
%                            of the gradient of u (0x0 is possible if the 
%                            gradient of u needs no further parameters)
%         parAlpha         - 'double' parameter alpha necessary for the 
%                            computation of the energy
%         minNrDof         - 'uint64' minimal number of dofs of the finest, red 
%                            refined, mesh on which the energy is approximated
%         minPrecision     - 'uint64' minimal number of significant digits the 
%                            approximation of the energy should possess
%         degree4Integrate - 'uint64' up to which the integration in integrate
%                            must be exact
%
% output: -

  addpath(genpath(pwd), genpath('../utils/'));

  if nargin < 11
    degree4Integrate = 20;
    if nargin < 10
      minPrecision = 2;
      if nargin < 9
        minNrDof = 1e4;
      end
    end
  end

  f = @(x) feval(fStr, x, fStrParams); 
  u = @(x) feval(uStr, x, uStrParams); 
  gradU = @(x) feval(gradUStr ,x, gradUStrParams);

  if strcmp(geometry, 'Polygon')
    [c4n, n4e, n4sDb, n4sNb] = ...
      computeGeometryPolygon(0);
    polygonMesh = true;
  else
    [c4n, n4e, n4sDb, n4sNb] = ...
      loadGeometry(geometry, 0);
    polygonMesh = false;
  end

  rhsStr = sprintf('alpha_%.30g_with_rhs_%s%s', ...
    parAlpha, fStr, sprintf('_%.30g', fStrParams));
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf(...
    'knownExactEnergies/%s/%s', ...
    geometry, rhsStr);
  mkdir(dirName);
  warning('on', 'MATLAB:MKDIR:DirectoryExists');

  output = struct;

  nrDof = [];
  energy = [];
  significantDigits = [0];
  
  while true
    % compute geometry
    [c4n, n4e, n4sDb, n4sNb] = refineUniformRed(c4n, n4e, n4sDb, n4sNb);
    if polygonMesh
      temp = unique(n4sDb);
      c4n(temp, :) = ...
        c4n(temp, :)./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
    end

    area4e  = computeArea4e(c4n, n4e);
    
    % compute nrDof
    s4e = computeS4e(n4e);
    n4s = computeN4s(n4e);
    tempStruct.n4sDb = n4sDb;
    tempStruct.s4n = computeS4n(n4e, n4s);
    tempStruct.nrSides = max(max(s4e));
    dof = computeDofCR(tempStruct);
    nrDof(end+1,1) = length(dof);
    output.nrDof = nrDof;

    % compute energy
    energy(end+1,1) = sum(...
        integrate(@(n4p, Gpts4p, Gpts4ref)(...
      parAlpha/2*u(Gpts4p).^2 + sqrt(sum(gradU(Gpts4p).^2, 2)) ...
      - f(Gpts4p).*u(Gpts4p)), ...
      c4n, n4e, degree4Integrate+1, area4e));
    output.energy = energy;

    % compute significant digits
    if length(energy) > 1 && fix(energy(end-1)) == fix(energy(end))
      % 20 is too precise, so it's sufficient
      dec1 = extractAfter(num2str(energy(end-1), '%.20f'), '.');
      dec2 = extractAfter(num2str(energy(end), '%.20f'), '.');
      for j = 1:min(length(dec1), length(dec2))
        if ~strcmp(dec1(j), dec2(j))
          break
        end
      end
      % in MATLAB j is known even after the loop ends
      significantDigits(end+1,1) = j-1;
    end
    output.significantDigits = significantDigits;

    % display status of computation
    fprintf([repmat('\n',1,50), ...
      '    minNrDof = %e\n    minPrecision = %d\n\n'], minNrDof, minPrecision);
    disp(struct2table(output));

    % save results
    name = sprintf('%s/nrDof_%d_significantDigits_%d.txt', ...
      dirName, nrDof(end), significantDigits(end));

    file = fopen(name, 'w');
    fprintf(file, 'nrDof   energy   significantDigits\n');
    fprintf(file, '%d   %.30g   %d\n', [nrDof, energy, significantDigits]');
    fclose(file);

    % check termination
    if nrDof(end) > minNrDof && significantDigits(end)>=minPrecision
      break
    end
  end

  % TODO maybe delete all but the best file with the name
  %  sprintf('%s/minPrecision_%d_nrDof_%d.txt', ...
  %    dirName, significantDigits(end), nrDof(end));
  %
  %    but this is really unimportant and probably not worth the effort


  % nrElems = size(n4e,1);
  % nrSides = max(max(s4e));
  % [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);
  % temp = zeros(nrSides,1);
  % for elem = 1 : nrElems
  %   temp(s4e(elem,:)) = ...
  %     temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
  % end
  % [~,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);

  % du = computeGradientNC(c4n,n4e,u);
  % E = computeEnergy(area4e,u,du,alpha,temp,MAMANC); end
