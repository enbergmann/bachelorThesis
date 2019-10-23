function exactEnergy(geometry, fStr, fStrParams, uStr, uStrParams, ...
    gradUStr, gradUStrParams, minNrDof, degree4Integrate, parAlpha, parBeta)
  % TODO write interface documentation
  % compute all energies until at least minNrDof
% Computes the degrees of freedom for the CR_0^1 FEM of the triangulation given
% by [n4e, n4sDb, n4sNb].
%
% computeDofCR.m
% input:  currData - struct with fields:
%                        n4sDb: nodes for Dirichlet boundary sides
%                          s4n: sides for nodes
%                      nrSides: number of sides
%
% output: dof      - '(1 x nrDof)-dimensional double array' where the j-th row 
%                    contains the number of the j-th degree of freedom


% TODO mode parameter? either compute until at least minNrDofs Dofs or
% until at least significantDigits significant digits (change from one mesh 
% to next smaller than 1e-'significant digits' or sth)

  %TODO func2str useless, probably just build string myself (%[rhs]_%[param1]_..._%paramN)

  addpath(genpath('../'), genpath('../../../utils/'));

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

  %TODO build str rhsStr here

  warning('off', 'MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf(...
    'knownExactEnergies/%s/%s', ...
    geometry, rhsStr);
  mkdir(dirName);
  warning('on',' MATLAB:MKDIR:DirectoryExists');
  name = sprintf('%s/%s.txt', minNrDof); %TODO  Write highest used nrDofs
                                          % specifically
        %TODO that means save everything to nrDofVec and energyVec and save 
        %those in the end (to be able to write highest nrDofs in the title)
  file = fopen(name, 'w');
  fprintf(file, 'nDoF   energy\n');

  for level = 0:red
    [c4n, n4e, n4sDb, n4sNb] = refineUniformRed(c4n, n4e, n4sDb, n4sNb);
    if polygonMesh
      temp = unique(n4sDb);
      c4n(temp, :) = ...
        c4n(temp, :)./repmat(sqrt(c4n(temp, 1).^2 + c4n(temp, 2).^2), 1, 2);
    end

    area4e = computeArea4e(c4n,n4e);
    
    energy = sum(...
        integrate(@(n4p,Gpts4p,Gpts4ref)(...
      alpha/2*u(Gpts4p).^2 + sqrt(sum(gradU(Gpts4p).^2,2)) - f(Gpts4p).*u(Gpts4p))...
      ,c4n,n4e,degree+1,area4e))
    

    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));
    dof = computeDof(n4e,nrSides,n4sDb,n4sNb);
    nrDof = length(dof);

    fprintf(file, '%d   %.30g\n',nrDof,energy);
  end

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
  % E = computeEnergy(area4e,u,du,alpha,temp,MAMANC);
end
