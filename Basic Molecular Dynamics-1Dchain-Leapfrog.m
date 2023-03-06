function [t u]=Leapfrog(ffunc,dt,N,y0)
%Leapfrog algorithm for the one-dimensional many-particle
%Hamiltonian dynamics
%Input variables:
% ffunc - force function
% dt - time step
% N - number of time steps
% y0 - a row vector of initial positions and momenta: y0 = [q0 p0]
%Output variables:
% t - a column vector of time points
% u - an N+1 by 2*M matrix of positions and momenta

t(1,1)=0; %t_0=0
M=length(y0)/2;
q=y0(1:M); %initial positions
p=y0(M+1:2*M); %initial momenta

for n=1:N,
    t(n+1,1)=t(n,1)+dt; %new time point: t_n+dt
    a=ffunc(q(n,:)); %forces at t_n
    p(n+1,:)=p(n,:)+dt*a; %momenta at t_n+dt/2
    q(n+1,:)=q(n,:)+p(n+1,:)*dt; %positions at t_n+dt
end

u=[q p];

return;