function plotAxisNC(c4n,n4e,u)
% TODO change to plotAxisCR
n4s = computeN4s(n4e);
nrSides = size(n4s,1);
mid4s = computeMid4s(c4n, n4s);

temp1=zeros(nrSides,1);
temp2=temp1;

temp1(mid4s(:,1)==0)=1;
temp2(mid4s(:,2)==0)=1;

[y,I1]=sort(mid4s(temp1==1,2));
[x,I2]=sort(mid4s(temp2==1,1));

u1=u(temp1==1);
u2=u(temp2==1);

u1=u1(I1);
u2=u2(I2);

plot(y,u1,'-*')
hold on
plot(x,u2,'-d')
