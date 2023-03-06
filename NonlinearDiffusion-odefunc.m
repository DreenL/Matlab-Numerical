function f=odefunc(t,y,h)

N=sqrt(length(y));
Y=reshape(y,[N N]);

F=zeros(N,N);
for i=1:N,
    ip=i+1-(i==N)*N;
    im=i-1+(i==1)*N;
    for j=1:N,
        jp=j+1-(j==N)*N;
        jm=j-1+(j==1)*N;
        uc=Y(i,j);
        ue=Y(ip,j);
        uw=Y(im,j);
        un=Y(i,jp);
        us=Y(i,jm);  
        te=mobfunc((uc+ue)/2)*(ue-uc);
        tw=mobfunc((uc+uw)/2)*(uw-uc);
        tn=mobfunc((uc+un)/2)*(un-uc);       
        ts=mobfunc((uc+un)/2)*(us-uc);        
        F(i,j)=(te+tw+tn+ts)/h^2;
    end
end

f=reshape(F,[N^2 1]);

return;

function M=mobfunc(phi)

%M=1;
M=(phi*(1-phi))^2;

return;