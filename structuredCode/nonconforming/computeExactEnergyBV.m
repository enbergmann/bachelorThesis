function computeExactEnergyBV(geometry, fStr, fStrParams, uStr, uStrParams, ...
    gradUStr, gradUStrParams, parAlpha, ...,
    minNrInnerEdges, minPrecision, degree4Integrate)
%% DOC
% Computes and saves an approximation (at least up to precision minPrecision
% and minNrInnerEdges inner edges) of the exact energy of a H^1_0 function u,
% whose pointwise gradient is known, for an input signal f on a mesh given by
% geometry.
%
% computeExactEnergyBV.m
% input: geometry         - 'char array with exactly one row' containing the 
%                           name of the geometry the user wants to approximate
%                           the exact energy on
%        fStr             - 'string'/'char array with exactly one row' 
%                           containing the name of the input signal f
%        fStrParams       - 'double array with exactly one row' containing the 
%                           necessary parameters to produce a function handle
%                           of f (0x0 is possible if f needs no further
%                           parameters)
%        uStr             - 'string'/'char array with exactly one row' 
%                           containing the name of the function u whose exact
%                           BV energy is to be approximated 
%        uStrParams       - 'double array with exactly one row' containing the 
%                           necessary parameters to produce a function handle
%                           of u (0x0 is possible if u needs no further
%                           parameters)
%        gradUStr         - 'string'/'char array with exactly one row' 
%                           containing the name of the gradient of the
%                           function u whose exact BV energy is to be
%                           approximated
%        gradUStrParams   - 'double array with exactly one row' containing the 
%                           necessary parameters to produce a function handle 
%                           of the gradient of u (0x0 is possible if the 
%                           gradient of u needs no further parameters)
%        parAlpha         - 'double' parameter alpha necessary for the 
%                           computation of the energy
%        minNrInnerEdges  - 'uint64' minimal number of inner edges of the
%                           finest, red refined, mesh on which the energy is
%                           approximated (optional)
%        minPrecision     - 'uint64' minimal number of significant digits the 
%                           approximation of the energy should possess
%                           (optional)
%        degree4Integrate - 'uint64' up to which the integration in integrate
%                           must be exact (optional)

%% INIT
  addpath(genpath(pwd), genpath('../utils/'));

  if nargin < 11
    degree4Integrate = 20;
    if nargin < 10
      minPrecision = 2;
      if nargin < 9, minNrInnerEdges = 1e4; end
    end
  end

  if isempty(fStrParams), f = @(x) feval(fStr, x); 
  else, f = @(x) feval(fStr, x, fStrParams); end
  if isempty(uStrParams), u = @(x) feval(uStr, x); 
  else, u = @(x) feval(uStr, x, uStrParams); end
  if isempty(gradUStrParams), gradU = @(x) feval(gradUStr, x);
  else, gradU = @(x) feval(gradUStr, x, gradUStrParams); end

  if strcmp(geometry, 'Polygon')
    [c4n, n4e, n4sDb, n4sNb] = computeGeometryPolygon(0);
    polygonMesh = true;
  else
    [c4n, n4e, n4sDb, n4sNb] = loadGeometry(geometry, 0);
    polygonMesh = false;
  end

  inSiStr = sprintf('alpha_%.30g_with_inSi_%s%s', ...
    parAlpha, fStr, sprintf('_%.30g', fStrParams));
  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf('knownExactEnergies/%s/%s', geometry, inSiStr);
  mkdir(dirName);
  warning('on', 'MATLAB:MKDIR:DirectoryExists');

  output = struct;

  nrInnerEdges = [];
  energy = [];
  significantDigits = 0;
  
%% MAIN
  while true
    % compute geometry
    [c4n, n4e, n4sDb, n4sNb] = refineUniformRed(c4n, n4e, n4sDb, n4sNb);
    if polygonMesh
      temp = unique(n4sDb);
      c4n(temp, :) = ...
        c4n(temp, :)./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
    end

    area4e  = computeArea4e(c4n, n4e);
    
    % compute nrInnerEdges
    s4e = computeS4e(n4e);
    n4s = computeN4s(n4e);
    tempStruct.n4sDb = n4sDb;
    tempStruct.s4n = computeS4n(n4e, n4s);
    tempStruct.nrSides = max(max(s4e));
    dof = computeDofCR(tempStruct);
    nrInnerEdges(end + 1, 1) = length(dof); %#ok<AGROW>
    output.nrInnerEdges = nrInnerEdges;

    % compute energy
    energy(end + 1, 1) = sum(...
        integrate(c4n, n4e, @(n4p, Gpts4p, Gpts4ref)(...
      parAlpha/2*u(Gpts4p).^2 + sqrt(sum(gradU(Gpts4p).^2, 2)) ...
      - f(Gpts4p).*u(Gpts4p)), degree4Integrate + 1, area4e)); %#ok<AGROW>
    output.energy = energy;

    % compute significant digits
    if length(energy) > 1 
      if fix(energy(end - 1))==fix(energy(end))
        dec1 = extractAfter(num2str(energy(end - 1), '%.20f'), '.');
        dec2 = extractAfter(num2str(energy(end), '%.20f'), '.');
          % 20 is too precise, so it's sufficient
        for j = 1:min(length(dec1), length(dec2))
          if ~strcmp(dec1(j), dec2(j)), break; end
        end
        significantDigits(end + 1, 1) = j - 1; %#ok<AGROW>
          % in MATLAB j is known even after the loop ends
      else
        significantDigits(end + 1, 1) = 0; %#ok<AGROW>
      end
    end
    output.significantDigits = significantDigits;

    % display status of computation
    fprintf([repmat('\n', 1, 50), ...
      '    minNrInnerEdges = %e\n    minPrecision = %d\n\n'], ...
      minNrInnerEdges, minPrecision);
    disp(struct2table(output));

    % save results
    name = sprintf('%s/nrInnerEdges%d_significantDigits_%d.txt', ...
      dirName, nrInnerEdges(end), significantDigits(end));

    file = fopen(name, 'w');
    fprintf(file, 'nrInnerEdges   energy   significantDigits\n');
    fprintf(file, '%d   %.30g   %d\n', ...
      [nrInnerEdges, energy, significantDigits]');
    fclose(file);

    % check termination
    if nrInnerEdges(end) > minNrInnerEdges && ...
        significantDigits(end) >= minPrecision, break; end
  end
end
