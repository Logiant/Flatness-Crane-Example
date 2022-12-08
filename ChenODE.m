function dydt = ChenODE(t, y, data)

%system parameters
m = data.m;
M = data.M;
g = data.g;
l = data.l;

F = data.u(t);

%state is [position, angle] --- [p, pDot, th, thDot]
pDot = y(2); th = y(3); thDot = y(4);

%explicitly calculate acceleration from the ODE
pDdot = (F + m*l*thDot^2*sin(th) + m*g*sin(th)) / (m + M - m*cos(th)); %control input

%use the rotational equation that is a function of pDot
thDdot = -(m*g*l*sin(th) + m*l*cos(th)*pDdot) / (m*l^2); %swing dynamics


dydt = [pDot; pDdot; thDot; thDdot];
end

