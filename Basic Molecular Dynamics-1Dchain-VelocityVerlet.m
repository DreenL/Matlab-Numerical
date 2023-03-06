function [t u]=VelocityVerlet(ffunc,dt,N,y0)
%Velocity Verlet algorithm for the one-dimensional many-particle
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
M=length(y0)/2; %number of particles
q=y0(1:M); %initial condition for the positions
p=y0(M+1:2*M); %intial condition for the momenta

a=ffunc(q); %forces at t=0
for n=1:N,
    t(n+1,1)=t(n,1)+dt; %next time point
    q(n+1,:)=q(n,:)+p(n,:)*dt+a*dt^2/2; %positions at t_n+dt
    an=ffunc(q(n+1,:)); %forces at t_n+dt
    p(n+1,:)=p(n,:)+(a+an)/2*dt; %momenta at t_n+dt
    a=an; %storing forces at t_n+dt for the next time step to avoid double calculation
end

u=[q p];

return;