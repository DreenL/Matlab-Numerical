function [t u]=RungeKutta(f,t0,T,y0,N,method)

 switch method
     case "ForwardEuler"
         s=1;
         c=[0];
         b=[1];
         A=[0];
     case "ExplicitMidpoint"
         s=2;
         c=[0 1/2];
         b=[0 1];
         A=[0 0; 1/2 0];
     case "Heun"
         s=2;
         c=[0 1];
         b=[1/2 1/2];
         A=[0 0; 1 0];
     case "Kutta3"
         s=3;
         c=[0 1/2 1];
         b=[1/6 2/3 1/6];
         A=[0 0 0; 1/2 0 0; -1 2 0];
     case "SSPRK3"
         s=3;
         c=[0 1 1/2];
         b=[1/6 1/6 2/3];
         A=[0 0 0; 1 0 0; 1/4 1/4 0];
 end
 
 [t u]=RungeKuttaGeneral(f,t0,T,y0,N,s,c,b,A);

 return;
         
         
    