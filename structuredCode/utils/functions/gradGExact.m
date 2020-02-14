function val = gradGExact(x, params)

parAlpha = params(1, 1);
parBeta = params(1, 2);

nP = size(x, 1);
val = zeros(nP, 2);
r = sqrt(sum(x.^2, 2));

temp = zeros(nP, 1);
% TODO not definied in r_temp = 0, hence val remains 0 there for now
% but choose a better value in [0,0] analyticaly
temp(0<r & r<=1/6) = 1;
%temp(r<=1/6) = 1;
r_temp = r(temp==1);
val(temp==1, :) = 108./r_temp.*x(temp==1, :);

temp = zeros(nP, 1);
temp(1/6<r & r<=1/3) = 1;
r_temp = r(temp==1);
val(temp==1, :) = (6*parAlpha*parBeta*(6*r_temp-1).^(parBeta-1) + ...
  1./r_temp.^3)./r_temp.*x(temp==1, :);

temp = zeros(nP, 1);
temp(1/3<r & r<=1/2) = 1;
r_temp = r(temp==1);
val(temp==1, :) = (cos(pi*(6*r_temp-2)).*(36*pi^2 + 1./r_temp.^2) + ...
  6*pi./r_temp.*sin(pi*(6*r_temp-2)))./r_temp.*x(temp==1, :);

temp = zeros(nP, 1);
temp(1/2<r & r<=5/6) = 1;
r_temp = r(temp==1);
val(temp==1, :) = (-6*parAlpha*parBeta*(5/2-3*r_temp).^(parBeta-1) - ...
  1./r_temp.^2)./r_temp.*x(temp==1, :);

temp = zeros(nP, 1);
temp(5/6<r & r<=1) = 1;
r_temp = r(temp==1);
val(temp==1, :) = -(18*pi^2*cos(pi*(6*r_temp-5)) + ...
  1./(4*r_temp.^2).*(1+cos(pi*(6*r_temp-5))) + ...
  3*pi./r_temp.*sin(pi*(6*r_temp-5)))./r_temp.*x(temp==1, :);

%val(sqrt(sum(x.^2,2))<.2) = 1; % abs(x)<0.2
