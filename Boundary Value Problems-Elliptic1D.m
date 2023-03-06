function [x u]=Elliptic1D(rhsfunc,range,N,bcmat)

a=range(1);
b=range(2);
h=(b-a)/N;
x=(a+h/2:h:b-h/2)';

alpha=bcmat(1,1);
beta=bcmat(1,2);
gamma=bcmat(1,3);

Aa=-(alpha/2+beta/h)/(alpha/2-beta/h);
Ba=-gamma/(alpha/2-beta/h);

alpha=bcmat(2,1);
beta=bcmat(2,2);
gamma=bcmat(2,3);

Ab=-(alpha/2-beta/h)/(alpha/2+beta/h);
Bb=-gamma/(alpha/2+beta/h);

L=(-2*diag(ones(1,N),0)+diag(ones(1,N-1),-1)+diag(ones(1,N-1),+1));
f=h^2*rhsfunc(x);
b=zeros(N,1);

L(1,1)=L(1,1)+Aa;
L(N,N)=L(N,N)+Ab;

b(1)=Ba;
b(N)=Bb;

u=L\(f+b);

return;