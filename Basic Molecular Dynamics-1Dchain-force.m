function a=force(q)
%force calculator for the one-dimensional particle chain
%Input variables:
% q - a row vector of particle positions
%Output variables:
% a - a row vector of forces

x=q'-q; %position difference matrix: x(i,j)=q(i)-q(j)
f=pforce(x); %force matrix: f(i,j)=f(q_i-q_j)
m=length(x); %number of particles
pos=(m*(0:m-1))+(1:m); %linear indices of the diagonal elements of 'f'
f(pos)=0; %setting the diagonal elements of 'f' to zero to avoid self-interactions
a=sum(f'); %the resulting force vector is a row vector: a(i) = sum_{j=1}^m f(i,j)

return;

function f=pforce(x)
f=12*x.^(-7).*(x.^(-6)-1);   
return;
    
