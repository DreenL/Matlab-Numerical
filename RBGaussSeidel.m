function [x y u err]=RBGaussSeidel_basic(rhsfunc,L,N,u0,numit)

h=L/N;
s=h/2:h:L-h/2;
[x y]=meshgrid(s,s);
x=x'; y=y';
f=rhsfunc(x,y);

u=u0;

for n=1:numit
    for i=1:N,
        ip=i+1-(i==N)*N;
        im=i-1+(i==1)*N;            
        for j=1:N,
            jp=j+1-(j==N)*N;
            jm=j-1+(j==1)*N;
            if ((mod(n,2)==0)&(mod(i+j,2)==0))|((mod(n,2)==1)&(mod(i+j,2)==1)),
                u(i,j)=(u(ip,j)+u(im,j)+u(i,jp)+u(i,jm)-h^2*f(i,j))/4;
            end
        end
    end
    tmp=0;
    for i=1:N,
        ip=i+1-(i==N)*N;
        im=i-1+(i==1)*N;            
        for j=1:N,
            jp=j+1-(j==N)*N;
            jm=j-1+(j==1)*N;  
            r=(u(ip,j)+u(im,j)+u(i,jp)+u(i,jm)-4*u(i,j))/h^2-f(i,j);
            tmp=tmp+r*r;
        end
    end
    err(n)=sqrt(tmp/N^2);
end

x=x'; y=y'; 
u=u';

return;