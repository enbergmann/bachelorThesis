function val = gradGExact(x, params)
  % TODO and for other functions: here x/|x| can be done after everything else
  % since it's in all definitions. should be more efficient
  % TODO still need to change some names now
  %
  % NOTE most of the time the pol variant is about 0.004s faster (by timeit())
  % TODO test this for f and uExact, since those are called very oftnoen
  % BUT profile with 500 calls of each function shows a difference of nearly
  % 3 secondes in favor of gradGExact, so this might be the better function
  % right here
  parAlpha = params(1, 1);
  parBeta = params(1, 2);
  
  r = sqrt(sum(x.^2, 2));
  val = zeros(length(r), 1);
  
  val(0<r & r<=1/6) = 108;
  
  ind = 1/6<r & r<=1/3;
  rTemp = r(ind);
  val(ind) = 6*parAlpha*parBeta*(6*rTemp-1).^(parBeta-1) + 1./rTemp.^2;
  
  ind = 1/3<r & r<=1/2;
  rTemp = r(ind);
  val(ind) = cos(pi*(6*rTemp-2)).*(36*pi^2 + 1./rTemp.^2) + ...
    6*pi./rTemp.*sin(pi*(6*rTemp-2));

  ind = 1/2<r & r<=5/6;
  rTemp = r(ind);
  val(ind) = -6*parAlpha*parBeta*(5/2-3*rTemp).^(parBeta-1) -  1./rTemp.^2;

  ind = 5/6<r & r<=1;
  rTemp = r(ind);
  val(ind) = -(18*pi^2*cos(pi*(6*rTemp-5)) + ...
    1./(2*rTemp.^2).*(1+cos(pi*(6*rTemp-5))) + ...
    3*pi./rTemp.*sin(pi*(6*rTemp-5)));
  
  % TODO not definied in rTemp = 0, hence val remains 0 there for now
  % but choose a better value in [0,0] analyticaly
  valTemp = zeros(length(r), 2);
  ind = r>0;
  rTemp = r(ind);
  valTemp(ind, :) = val(ind)./rTemp.*x(ind, :);
  val = valTemp;
end
