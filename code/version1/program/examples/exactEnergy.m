function exactEnergy(red,degree,alpha,delta)
  addpath(genpath('../'),genpath('../../../utils/'));
  f = @(x) g(x,alpha,delta);  
  u = @(x) gUexact(x,alpha,delta);
  gradU = @(x) gradGuExact(x,delta);
  % [c4n,n4e] = computeGeometryPolygon(red);
  %[c4n,n4e] = computeGeometryCriss(red);
  [c4n,n4e] = computeGeometryBigSquare(red);
  area4e = computeArea4e(c4n,n4e);
  
  energy = sum(...
      integrate(@(n4p,Gpts4p,Gpts4ref)(...
    alpha/2*u(Gpts4p).^2 + sqrt(sum(gradU(Gpts4p).^2,2)) - f(Gpts4p).*u(Gpts4p))...
    ,c4n,n4e,degree+1,area4e))
  
   
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
