function f=odefunc(t,y,h,alpha)

N=sqrt(length(y));
Y=reshape(y,[N N]);

beta=(1+sqrt(5))/4*alpha;

b1=(1-beta)^2;
b2=beta*(1-beta);
b3=beta^2;

L=[1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1];
s1=2+1/alpha;
s2=2-1/alpha;
K1=[s1 0 s2 0; 0 s1 0 s2; s2 0 s1 0; 0 s2 0 s1];
C=(sqrt(5)-1)/2;
h1=(1+C/alpha)/2;
h2=(1-C/alpha)/2;
K2=[h1 0 h2 0; 0 h1 0 h2; h2 0 h1 0; 0 h2 0 h1];

F=zeros(N,N);
for i=1:N,
    ip=i+1-(i==N)*N;
    im=i-1+(i==1)*N;
    for j=1:N,
        jp=j+1-(j==N)*N;
        jm=j-1+(j==1)*N;
        v=[Y(ip,j) Y(i,jp) Y(im,j) Y(i,jm)]';
        phi(1:4,:)=(1-alpha)*Y(i,j)+alpha*v;
        w=[Y(ip,jp) Y(im,jp) Y(im,jm) Y(ip,jm)]';        
        phi(5:8,:)=b1*Y(i,j)+b2*L*v+b3*w;
        M=mobfunc(phi);
        c1=K1*M(1:4);
        c2=K2*M(5:8);
        F(i,j)=(dot(c1,v)+dot(c2,w)-(4*sum(M(1:4))+sum(M(5:8)))*Y(i,j))/(6*h^2);
    end
end

f=reshape(F,[N^2 1]);

return;

function M=mobfunc(phi)

M=(phi.*(1-phi)).^2;

return;
