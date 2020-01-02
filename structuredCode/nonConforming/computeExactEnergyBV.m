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
%         

%TODO (maybe)
%save every step in file, but then it's hard to fix the problem that the names
%of incomplete runs might coincide with complete, and therefore, better runs
%
%for now this is fine, but it needs to complete the computation to yield a 
%result, which isn't optimale but the servers might be able to handle it

  addpath(genpath(pwd), genpath('../utils/'));

  if nargin < 12
    degree4Integrate = 20;
    if nargin < 11
      minPrecision = 3;
      if nargin < 10
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

  rhsStr = sprintf('%s%s', fStr, sprintf('_%.30g', fStrParams));
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf(...
    'knownExactEnergies/%s/%s', ...
    geometry, rhsStr);
  mkdir(dirName);
  warning('on', 'MATLAB:MKDIR:DirectoryExists');

  nrDofVec = [];
  energyVec = [0];

  while true
    [c4n, n4e, n4sDb, n4sNb] = refineUniformRed(c4n, n4e, n4sDb, n4sNb);
    if polygonMesh
      temp = unique(n4sDb);
      c4n(temp, :) = ...
        c4n(temp, :)./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
    end

    area4e  = computeArea4e(c4n, n4e);
    
    energyVec(end+1) = sum(...
        integrate(@(n4p, Gpts4p, Gpts4ref)(...
      parAlpha/2*u(Gpts4p).^2 + sqrt(sum(gradU(Gpts4p).^2, 2)) ...
      - f(Gpts4p).*u(Gpts4p)), ...
      c4n, n4e, degree4Integrate+1, area4e));
    
    s4e = computeS4e(n4e);
    n4s = computeN4s(n4e);
    tempStruct.n4sDb = n4sDb;
    tempStruct.s4n = computeS4n(n4e, n4s);
    tempStruct.nrSides = max(max(s4e));
    dof = computeDofCR(tempStruct);
    nrDofVec(end+1) = length(dof);
    if nrDofVec(end) > minNrDof ...
        && abs(energyVec(end)-energyVec(end-1)) < 10^(-minPrecision)
      break
    end

    fprintf('---------------------------------------\nenergy  nrDof\n');
    fprintf('%.20g  %d\n', energyVec(end), nrDofVec(end));
  end

  name = sprintf('%s/minPrecision_%d_nrDof_%d.txt', dirName, minPrecision, nrDofVec(end));

  file = fopen(name, 'w');
  fprintf(file, 'nrDof   energy\n');
  fprintf(file, '%d   %.30g\n', [nrDofVec; energyVec(2:end)]);
  fclose(file);
   


  % nrElems = size(n4e,1);
  % nrSides = max(max(s4e));
  % [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);
  % temp = zeros(nrSides,1);
  % for elem = 1 : nrElems
  %   temp(s4e(elem,:)) = temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
  % end
  % [~,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);

  % du = computeGradientNC(c4n,n4e,u);
  % E = computeEnergy(area4e,u,du,alpha,temp,MAMANC); end
