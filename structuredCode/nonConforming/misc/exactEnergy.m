function exactEnergy(geometry, rightHandSide, paramsRHS, exactSolution, paramsSol, ...
    minNrDof, degree4Integrate, parAlpha, parBeta)
  % TODO
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
  % compute all energies until at least minNrDof

  %TODO idea feval function names with strArray paramsList and later
  %use func2str to create a folder for the results
  %TODO splice param list and insert params in feval (and think about what to do
  %when it's empty)
  %TODO test if the naming with func2str yields sth like @(x)g(x,1,1)
  addpath(genpath('../'),genpath('../../../utils/'));
  %f = @(x) g(x,alpha,delta);  
  %u = @(x) gUexact(x,alpha,delta);
  %gradU = @(x) gradGuExact(x,delta);
  f = feval(rightHandSide, paramsRHS) %TODO this is just a sketch

  % TODO cut

  [c4n,n4e,n4sDb,n4sNb] = computeGeometryPolygon(0);
  % [c4n,n4e] = computeGeometryCriss(red);
  % [c4n,n4e,n4sDb,n4sNb] = computeGeometryBigSquare(0);

  %dirName = sprintf(...
  %  '../../../../results/exactEnergy/nonconforming/%s/%s',...
  %  geometry,datestr(now,'yy_mm_dd_HH_MM_SS'));
  warning('off','MATLAB:MKDIR:DirectoryExists');
  dirName = sprintf(...
    'knownExactEnergies/%s/%s',...
    geometry,func2str(f));
  mkdir(dirName);
  warning('on','MATLAB:MKDIR:DirectoryExists');
  name = sprintf('%s/%s.txt',minNrDof); %TODO maybe write highest used nrDofs
                                          % specifically
  file = fopen(name,'w');
  fprintf(file, 'nDoF   energy\n');

  for level = 0:red
    [c4n,n4e,n4sDb,n4sNb] = refineUniformRed(c4n,n4e,n4sDb,n4sNb);
    if strcmp(geometry,'polygon')
      temp=unique(n4sDb);
      c4n(temp,:)=c4n(temp,:)./repmat(sqrt(c4n(temp,1).^2+c4n(temp,2).^2),1,2);
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
