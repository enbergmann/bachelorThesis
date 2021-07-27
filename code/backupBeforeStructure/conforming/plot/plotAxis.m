function plotAxis(c4n,u)
nC=size(c4n,1);

temp1=zeros(nC,1);
temp2=temp1;

temp1(c4n(:,1)==0)=1;
temp2(c4n(:,2)==0)=1;

[y,I1]=sort(c4n(temp1==1,2));
[x,I2]=sort(c4n(temp2==1,1));

u1=u(temp1==1);
u2=u(temp2==1);

u1=u1(I1);
u2=u2(I2);

plot(y,u1,'-*')
hold on
plot(x,u2,'-d')