function computeExactEnergyBV(geometry, fStr, fStrParams, uStr, uStrParams, ...
    gradUStr, gradUStrParams, parAlpha, parBeta, ...
    minNrDof, minPrecision, degree4Integrate)

% Computes the degrees of freedom for the CR_0^1 FEM of the triangulation given
% by [n4e, n4sDb, n4sNb].
%
% computeExactEnergyBV.m
% input:  geometry         - 'str'
%         fStr             -
%         fStrParams       -
%         uStr             -
%         uStrParams       -
%         gradUSr          -
%         gradUStrParams   -
%         parAlpha         -
%         parBeta          - 
%         minNrDof         - 
%         minPrecision     -
%         degree4Integrate -
%
% output: -
%         

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
