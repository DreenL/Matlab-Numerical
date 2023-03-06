function f=odefun(t,y)
%Right-hand side of the one dimensional M-particle Hamiltonian dynamics
%t is a scalar and y is a column vector
%The first M elements of y are the positions, and the second M elements are the
%momenta

M=length(y)/2; %number of particles

q=y(1:M); %positions in the first half of 'y'
p=y(M+1:2*M); %momenta in the second half of 'y'

%Setting up the first M elements of 'f' (column vector)
f(1:M,1)=p; %dq_i/dt = p_i (first half of 'f' = second half of 'y')

%force calculation for ALL particles
for i=1:M,
    dq=q(i)*ones(M,1)-q; %column vector of position differences: dq(j)=q(i)-q(j)
    dq(i)=[]; %dropping the 'i'th element of dq to avoid self-interaction
    F=pforce(dq); %column vector of forces on particle 'i' from other particles
    f(M+i,1)=sum(F); %total force on particle 'i'
end

return;

function f=pforce(x)
f=12*x.^(-7).*(x.^(-6)-1);   
return;
    
