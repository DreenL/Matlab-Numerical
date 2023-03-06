function [t u]=ForwardEuler(f,t0,T,y0,N)
%Forward Euler method to solve the scalar Cauchy problem
%Input variables:
%f - the right-hand side of the ODE f(t,y)
%t0 - initial time point
%T - length of the time interval: t0 <= t <= t0+T
%y0 - initial condition: y(t0)=y0
%N - number of sub-intervals (points) on the uniform grid
%Output variables:
%t - vector of grid points t0,t1,t2,...,tN
%u - vector of numerical solutions u0,u1,u2,...,uN 

dt=T/N; %time step
%Applying initial condition
t(1)=t0; %the first element of vector 't' is t0
u(1)=y0; %the first element of vector 'u' is y0

%time loop to calucalte u1,u2,...,uN
%MATLAB indexing runs from 1, therefore
%t(1)=t0, t(2)=t1, t(3)=t2,...,t(N+1)=tN, and
%u(1)=y0, u(2)=u1, u(3)=u2,...,u(N+1)=uN.
for n=1:N,
    tn=t(n); un=u(n); %temporary variables
    fn=f(tn,un); %value of f(t,y) at t=tn
    t(n+1)=tn+dt; %next grid point
    u(n+1)=un+dt*fn; %numerical solution at t(n+1)
end

%the output arrays 't' and 'u' have N+1 elements as:
%t = (t0,t1,t2,...,tN)
%u = (u0,u1,u2,...,uN)

return;