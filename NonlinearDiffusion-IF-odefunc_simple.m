function f=odefunc(t,y,h,alpha)

N=sqrt(length(y));
Y=reshape(y,[N N]);

beta=(1+sqrt(5))/4*alpha;

b1=(1-beta)^2;
b2=beta*(1-beta);
b3=beta^2;
C=(sqrt(5)-1)/2;

F=zeros(N,N);
for i=1:N,
    ip=i+1-(i==N)*N;
    im=i-1+(i==1)*N;
    for j=1:N,
        jp=j+1-(j==N)*N;
        jm=j-1+(j==1)*N;
        tmp=(1-alpha)*Y(i,j)+alpha*Y(ip,j); M(1)=mobfunc(tmp);
        tmp=(1-alpha)*Y(i,j)+alpha*Y(i,jp); M(2)=mobfunc(tmp);
        tmp=(1-alpha)*Y(i,j)+alpha*Y(im,j); M(3)=mobfunc(tmp);
        tmp=(1-alpha)*Y(i,j)+alpha*Y(i,jm); M(4)=mobfunc(tmp);
        tmp=b1*Y(i,j)+b2*(Y(ip,j)+Y(i,jp))+b3*Y(ip,jp); M(5)=mobfunc(tmp);
        tmp=b1*Y(i,j)+b2*(Y(im,j)+Y(i,jp))+b3*Y(im,jp); M(6)=mobfunc(tmp);        
        tmp=b1*Y(i,j)+b2*(Y(im,j)+Y(i,jm))+b3*Y(im,jm); M(7)=mobfunc(tmp);
        tmp=b1*Y(i,j)+b2*(Y(i,jm)+Y(ip,j))+b3*Y(ip,jm); M(8)=mobfunc(tmp);
        rhs=-(4*sum(M(1:4))+sum(M(5:8)))*Y(i,j);
        rhs=rhs+(2*(M(1)+M(3))+(M(1)-M(3))/alpha)*Y(ip,j);
        rhs=rhs+(2*(M(2)+M(4))+(M(2)-M(4))/alpha)*Y(i,jp);
        rhs=rhs+(2*(M(1)+M(3))+(M(3)-M(1))/alpha)*Y(im,j);
        rhs=rhs+(2*(M(2)+M(4))+(M(4)-M(2))/alpha)*Y(i,jm);
        rhs=rhs+(M(5)+M(7)+C*(M(5)-M(7))/alpha)/2*Y(ip,jp);
        rhs=rhs+(M(6)+M(8)+C*(M(6)-M(8))/alpha)/2*Y(im,jp);
        rhs=rhs+(M(5)+M(7)-C*(M(5)-M(7))/alpha)/2*Y(im,jm);
        rhs=rhs+(M(6)+M(8)-C*(M(6)-M(8))/alpha)/2*Y(ip,jm);        
        F(i,j)=rhs/(6*h^2);
    end
end

f=reshape(F,[N^2 1]);

return;

function M=mobfunc(phi)

M=(phi*(1-phi))^2;

return;
