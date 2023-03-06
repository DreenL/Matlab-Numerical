function E=energy(q,p)
%energy function for the one-dimensional M-particle system
% q - row vector of positions 
% p - row vector of momenta

K=dot(p,p)/2; %kinetic energy

%potential energy
x=q'-q; % position difference matrix: x(i,j)=q(i)-q(j)
u=potential(x); %potential energy on matrix x: u(i,j)=u(q_i-q_j)
m=length(q); %number of particles
pos=(m*(0:m-1))+(1:m); %linear indices of the diagonal elements of matrix 'u'
u(pos)=0; %setting the diagonal elements of 'u' to 0 to avoid self-interactions
U=sum(sum(u))/2; %adding all elements of 'u'

E=K+U;

return;

function u=potential(x)
u=x.^(-6).*(x.^(-6)-2);
return;