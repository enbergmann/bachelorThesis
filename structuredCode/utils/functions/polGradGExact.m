function val = polGradGExact(x, params)
  parAlpha = params(1, 1);
  parBeta = params(1, 2);
  
  [phi, r] = cart2pol(x(:, 1), x(:, 2));
  
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
  val(ind) = -6*parAlpha*parBeta*(5/2-3*rTemp).^(parBeta-1) - 1./rTemp.^2;

  ind = 5/6<r & r<=1;
  rTemp = r(ind);
  val(ind) = -(18*pi^2*cos(pi*(6*rTemp-5)) + ...
    1./(2*rTemp.^2).*(1+cos(pi*(6*rTemp-5))) + ...
    3*pi./rTemp.*sin(pi*(6*rTemp-5)));
  
  val = val.*[cos(phi), sin(phi)];
    % TODO phi not defined for x = [0, 0], find out what to do, for now ignore
    % those phi since val = 0 in those components as defined above
end
