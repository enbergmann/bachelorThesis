function val = gUexact(x,alpha,delta)

nP=size(x,1);
val = zeros(nP,1);
r=sqrt(sum(x.^2,2));

temp = zeros(nP,1);
temp(r<=1/6) = 1;
val(temp==1) = 1;

temp = zeros(nP,1);
temp(1/6<r & r<=1/3) = 1;
r_temp = r(temp==1);
val(temp==1) = 1 + (6*r_temp - 1).^delta;

temp = zeros(nP,1);
temp(1/3<r & r<=1/2) = 1;
val(temp==1) = 2;

temp = zeros(nP,1);
temp(1/2<r & r<=5/6) = 1;
r_temp = r(temp==1);
val(temp==1) = 2*(5/2 - 3*r_temp).^delta;

temp = zeros(nP,1);
temp(5/6<r & r<=1) = 1;
r_temp = r(temp==1);
val(temp==1) = 0;

%val(sqrt(sum(x.^2,2))<.2) = 1; % abs(x)<0.2
