function drawf(red)
  alpha = 1;
  delta = 1;
  [c4n,n4e,~,~] = computeGeometryPolygon(red); 
  u=g(c4n,alpha,delta);
  figure;
  trisurf(n4e,c4n(:,1),c4n(:,2),u);
  figure;
  plotAxis(c4n,u);
end
