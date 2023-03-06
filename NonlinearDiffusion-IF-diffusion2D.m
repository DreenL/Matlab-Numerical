clear all;
close all;

L=1;
N=64;
T=1/64;
alpha=0.75;

h=L/N;
s=h/2:h:L-h/2;
[x y]=meshgrid(s,s);

U0=exp(-32*((x-L/2).^2+(y-L/2).^2));

u0=reshape(U0,[N^2 1]);
%[t u]=ode23(@(t,y)odefunc_simple(t,y,h,alpha),[0 T],u0);
[t u]=ode23(@(t,y)odefunc_advanced(t,y,h,alpha),[0 T],u0);

nframes=6;
tstep=T/(nframes-1);
for i=1:nframes,
    tpos=(i-1)*tstep;
    [tmp idx]=min(abs(t-tpos));
    figure;
    contour(x,y,reshape(u(idx,:),[N N]),25);
    grid on;
    axis equal;
    surf(x,y,reshape(u(idx,:),[N N]));
    axis([0 L 0 L 0 max(max(U0))]);
    ftitle=sprintf("t = %f",t(idx));
    title(ftitle);
    xlabel('x');
    ylabel('y');
    zlabel('u');
    box on;
    set(gca,'fontsize',12);
end