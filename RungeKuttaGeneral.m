function [t u]=RungeKuttaGeneral(f,t0,T,y0,N,s,c,b,A)
%Runge-Kutta method to solve the scalar Cauchy problem
%Input variables:
% f - the right-hand side of the ODE f(t,y)
% t0 - initial time point
% T - length of the time interval: t0 <= t <= t0+T
% y0 - initial condition: y(t0)=y0
% N - number of sub-intervals (points) on the uniform grid
% s - number of stages
% c - 'c' vector of the Butcher array
% b - 'b' vector of the Butcher array
% A - 'A' matrix of the Butcher array
%Output variables:
%t - vector of grid points t0,t1,t2,...,tN
%u - vector of numerical solutions u0,u1,u2,...,uN 

dt=T/N; %time step
u(1)=y0; %first time point
t(1)=t0; %initial condition

K=zeros(1,s); %initialising the 'K' vector (dot(.,.) needs this!)

for n=1:N,
    %New time point
    t(n+1)=t(n)+dt;
    %Calculating the 'K_i' values 
    for i=1:s, %loop for the stages
        ts=t(n)+c(i)*dt; %'t' variable      
        us=u(n)+dt*dot(A(i,:),K); %'u' variable
        K(i)=f(ts,us); %calling the righ-hand side of the ODE
    end
    %Calculating u_{n+1} using the 'K' vector
    u(n+1)=u(n)+dt*dot(b,K);
end