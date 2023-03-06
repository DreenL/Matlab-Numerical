function [t u]=CrankNicolson(f,t0,T,y0,N,df)
%Crank-Nicolson method to solve the scalar Cauchy problem
%Input variables:
%f - the right-hand side of the ODE f(t,y)
%t0 - initial time point
%T - length of the time interval: t0 <= t <= t0+T
%y0 - initial condition: y(t0)=y0
%N - number of sub-intervals (points) on the uniform grid
%df - partial derivative of f(t,y) w.r.t. 'y'
%Output variables:
%t - vector of grid points t0,t1,t2,...,tN
%u - vector of numerical solutions u0,u1,u2,...,uN 

dt=T/N; %time step
%Applying initial condition
t(1)=t0; %the first element of vector 't' is t0
u(1)=y0; %the first element of vector 'u' is y0

%time loop to calucalte u1,u2,...,uN
for n=1:N,
    tn=t(n); un=u(n); %temporary variables
    fn=f(tn,un); %value of f(t,y) at t=tn
    t(n+1)=tn+dt; %next grid point
    z=un+dt*fn; %guess for u(n+1) using forward Euler method
    err=un+(dt/2)*(fn+f(t(n+1),z))-z; %error of the initial guess
    %Newton method to approximate u(n+1) more precisely  
    while abs(err)>1E-8,
        derr=(dt/2)*df(t(n+1),z)-1; %derivative of the Crank-Nicolson sheme rearranged to 0
        z=z-err/derr; %refining the numerical solution
        err=un+(dt/2)*(fn+f(t(n+1),z))-z; %error of the refined numerical solution
    end
    %As |err|<1E-8 at this point, we accept 'z' as u(n+1)
    u(n+1)=z;
end

return;