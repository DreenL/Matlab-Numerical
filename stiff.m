clear all;
close all;

t0=0;
T=300;
y0=[2 0];

[t1 u1]=ode45(@odefun,[t0 t0+T],y0);
figure;
plot(t1,u1(:,1),'o-');
grid on; box on;
xlabel('t');
ylabel('u^{(1)}');
title('ode45');
fprintf("ode45 needed %d steps.\n",length(t1));

[t2 u2]=ode15s(@odefun,[t0 t0+T],y0);
figure;
plot(t2,u2(:,1),'s-');
grid on; box on;
xlabel('t');
ylabel('u^{(1)}');
title('ode15s');
fprintf("ode15s needed %d steps.\n",length(t2));