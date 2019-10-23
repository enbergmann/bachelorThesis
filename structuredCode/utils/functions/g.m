function val = g(x, params)

alpha = params(1, 1);
delta = params(1, 2);

nP=size(x,1);
val = zeros(nP,1);
r=sqrt(sum(x.^2,2));

temp = zeros(nP,1);
temp(r<=1/6) = 1;
r_temp = r(temp==1);
val(temp==1) = alpha-12*(2-9*r_temp);

temp = zeros(nP,1);
temp(1/6<r & r<=1/3) = 1;
r_temp = r(temp==1);
val(temp==1) = alpha*(1+(6*r_temp-1).^delta)-1./r_temp;

temp = zeros(nP,1);
temp(1/3<r & r<=1/2) = 1;
r_temp = r(temp==1);
val(temp==1) = 2*alpha+6*pi*sin(pi*(6*r_temp-2))-1./r_temp.*cos(pi*(6*r_temp-2));

temp = zeros(nP,1);
temp(1/2<r & r<=5/6) = 1;
r_temp = r(temp==1);
val(temp==1) = 2*alpha*(5/2-3*r_temp).^delta+1./r_temp;

temp = zeros(nP,1);
temp(5/6<r & r<=1) = 1;
r_temp = r(temp==1);
val(temp==1) = -3*pi*sin(pi*(6*r_temp-5))+1./(2*r_temp).*(1+cos(pi*(6*r_temp-5)));

%val(sqrt(sum(x.^2,2))<.2) = 1; % abs(x)<0.2
