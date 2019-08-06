function E = exactDiscreteEnergy(c4n,n4e,s4e,u)
% Compute approximate exact discrete energy, u; coefficients of CR interpolation
  alpha = 1;
  delta = 1;
  f=@(x)g(x,alpha,delta);  
  nrElems = size(n4e,1);
  area4e = computeArea4e(c4n,n4e);
  nrSides = max(max(s4e));
  [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,200,area4e);
  temp = zeros(nrSides,1);
  for elem = 1 : nrElems
    temp(s4e(elem,:)) = temp(s4e(elem,:)) + [temp1(elem),temp2(elem),temp3(elem)]';
  end
  [~,MAMANC] = computeFeMatrices(c4n,n4e,s4e,area4e,nrElems);

  du = computeGradientNC(c4n,n4e,u);
  E = computeEnergy(area4e,u,du,alpha,temp,MAMANC);
end
