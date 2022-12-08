clear; close all; clc;

VIZ = false;

%system parameters
m = 200; M = 50; g = 9.81; l = 5;
alpha = 0.5; %weight for the control penalty term

%initial and final position
p0 = 0;
pf = 5;
%initial and final speed
v0 = 1.5;
vf = 0;
%initial and final theta
th0 = -pi/36;
thf = 0;
%initial and final angular speed
w0 = 0;
wf = 0;
%initial and final time
t0 = 0; 
tf = 15;

%constraints:
p_max = 10;
p_min = -10;

%solver parameters
nt = 1000;

%% concatenateinitial conditions
ICs = [p0 + l*th0; %yo
       v0 + l*w0; %dy0
       -g*th0; %ddy0
       -g*w0];%dddy0
FCs = [pf + l*thf; %yf
       vf + l*wf; %dyf
       -g*thf; %ddyf
       -g*wf];%dddyg

   
%define some useful shorthand variables
w = sqrt(g/l); %natural frequency of the pendulum
ag = alpha*g; %shows up in the exponentials
   
%build the parameter struct and ODE handle
data = struct('M', M, 'm', m, 'g', g, 'l', l, 'di', 1, ...
              'alpha', alpha, 'p_min', p_min, 'p_max', p_max);
   
% solve for the unconstrained trajectory
tic
C = UnconstrainedBVP(ICs, FCs, [t0, tf], data);
toc
%% generate the flat variable trajectories
t = linspace(t0, tf, nt);

y = C(1) + C(2)*t + C(3)*t.^2 + C(4)*t.^3 + C(5)*t.^4 + C(6)*t.^5 ...
    + C(7)*exp(t./(ag)) + C(8)*exp(-t./(ag));

dy = C(2) + 2*C(3)*t + 3*C(4)*t.^2 + 4*C(5)*t.^3 + 5*C(6)*t.^4 ...
    + C(7)*exp(t./(ag))/(ag) - C(8)*exp(-t./(ag))/(ag);

ddy = 2*C(3) + 6*C(4)*t + 12*C(5)*t.^2 + 20*C(6)*t.^3 ...
    + C(7)*exp(t./(ag))/(ag)^2 + C(8)*exp(-t./(ag))/(ag)^2;

dddy = 6*C(4) + 24*C(5)*t + 60*C(6)*t.^2 ...
    + C(7)*exp(t./(ag))/(ag)^3 - C(8)*exp(-t./(ag))/(ag)^3;

ddddy = 24*C(5) + 120*C(6)*t ...
    + C(7)*exp(t./(ag))/(ag)^4 + C(8)*exp(-t./(ag))/(ag)^4;

%% run the unconstrained simulation
F = (M+m)*(ddy + l/g*ddddy) + m*l*(ddy/g.*(dddy/g).^2 - ddddy/g);

th = -ddy/g;
dth = -dddy/g;

%generate the control action
u = griddedInterpolant(t, F);
data.u = u;

%state is [position, angle] --- [x, xDot, th, thDot]
x0 = [p0; v0; th0; w0];

%run and visualize the simulation
[t_sim, p_sim, v_sim, th_sim, w_sim] = ...
                        Simulation(x0, linspace(t0, tf, nt), data);
                    
                    
                    
                    
if VIZ
    Visualize(p_sim, th_sim, t_sim, data, 1);
end

lw = 2.75;
figure(12); clf; hold on;
plot([t0, tf], [p_max, p_max], '--r', 'linewidth', lw/2);
plot(t_sim, p_sim, 'linewidth', lw, 'color',  [0 0.4470 0.7410])


figure(11); clf; hold on;
plot(t_sim, th_sim, 'linewidth', lw, 'color',  [0 0.4470 0.7410])


%% Check what time the constraint becomes violated

p = y + l/g*ddy;
t1 = t(find(p > p_max, 1));

if isempty(t1)
    fprintf('No constraint violation occured!\n');
    return
end

fprintf('Constraint violation near t=%g s\n', t1);


%helper objective function to pass data to the optmizer
obj = @(x) Optimizer(x, t0, tf, ICs, FCs, data);


%initial guess
t_guess = [t1+1; t1+3];

%bound t1 and t2 between the initial and final times
lb = [t0; t0];
ub = [tf; tf];

%solve the optimization problem to get t1, t2
opts = optimoptions('fmincon', 'Display', 'None', ...
                     'Algorithm', 'interior-point'); %interior-point or sqp

tic
[x,J] = fmincon(obj, t_guess, [], [], [], [], lb, ub, @const, opts);
t1 = x(1); t2 = x(2);
toc
fprintf("Optimality gap of J = %g\n", J);

%check the BVP solution
C = ConstrainedBVP(t1, t2, t0, tf, ICs, FCs, data);

%add t1 and t2 explicitly
t = [t(t<t1), t1, t(t>t1&t<t2), t2, t(t>t2)];

%generate the first unconstrained arc
ta = t(t<=t1);

ya = C(1) + C(2)*ta + C(3)*ta.^2 + C(4)*ta.^3 + C(5)*ta.^4 + C(6)*ta.^5 ...
    + C(7)*exp(ta./(ag)) + C(8)*exp(-ta./(ag));

dya = C(2) + 2*C(3)*ta + 3*C(4)*ta.^2 + 4*C(5)*ta.^3 + 5*C(6)*ta.^4 ...
    + C(7)*exp(ta./(ag))/(ag) - C(8)*exp(-ta./(ag))/(ag);

ddya = 2*C(3) + 6*C(4)*ta + 12*C(5)*ta.^2 + 20*C(6)*ta.^3 ...
    + C(7)*exp(ta./(ag))/(ag)^2 + C(8)*exp(-ta./(ag))/(ag)^2;

dddya = 6*C(4) + 24*C(5)*ta + 60*C(6)*ta.^2 ...
    + C(7)*exp(ta./(ag))/(ag)^3 - C(8)*exp(-ta./(ag))/(ag)^3;

ddddya = 24*C(5) + 120*C(6)*ta ...
    + C(7)*exp(ta./(ag))/(ag)^4 + C(8)*exp(-ta./(ag))/(ag)^4;

%generate the constrained arc
tb = t(t>t1&t<t2);

yb = C(9)*cos(w*tb) + C(10)*sin(w*tb) + p_max;

dyb = -w*C(9)*sin(w*tb) + w*C(10)*cos(w*tb);

ddyb = -w^2*C(9)*cos(w*tb) - w^2*C(10)*sin(w*tb);

dddyb = w^3*C(9)*sin(w*tb) - w^3*C(10)*cos(w*tb);

ddddyb = w^4*C(9)*cos(w*tb) - w^4*C(10)*sin(w*tb);


%generate the second unconstrained arc
tc = t(t>=t2);

yc = C(11) + C(12)*tc + C(13)*tc.^2 + C(14)*tc.^3 + C(15)*tc.^4 + C(16)*tc.^5 ...
    + C(17)*exp(tc./(ag)) + C(18)*exp(-tc./(ag));

dyc = C(12) + 2*C(13)*tc + 3*C(14)*tc.^2 + 4*C(15)*tc.^3 + 5*C(16)*tc.^4 ...
    + C(17)*exp(tc./(ag))/(ag) - C(18)*exp(-tc./(ag))/(ag);

ddyc = 2*C(13) + 6*C(14)*tc + 12*C(15)*tc.^2 + 20*C(16)*tc.^3 ...
    + C(17)*exp(tc./(ag))/(ag)^2 + C(18)*exp(-tc./(ag))/(ag)^2;

dddyc = 6*C(14) + 24*C(15)*tc + 60*C(16)*tc.^2 ...
    + C(17)*exp(tc./(ag))/(ag)^3 - C(18)*exp(-tc./(ag))/(ag)^3;

ddddyc = 24*C(15) + 120*C(16)*tc ...
    + C(17)*exp(tc./(ag))/(ag)^4 + C(18)*exp(-tc./(ag))/(ag)^4;

y_c = [ya, yb, yc];
dy_c = [dya, dyb, dyc];
ddy_c = [ddya, ddyb, ddyc];
dddy_c = [dddya, dddyb, dddyc];
ddddy_c = [ddddya, ddddyb, ddddyc];

t_c = [ta, tb, tc];

%% calculate force, theta, etc

F_c = (M+m)*(ddy_c + l/g*ddddy_c) + m*l*(ddy_c/g.*(dddy_c/g).^2 - ddddy_c/g);

theta_c = -ddy_c./g;

figure(11); hold on;
plot(t_c, theta_c, '-.', 'linewidth', lw, 'color', [0.8500 0.3250 0.0980]);


p_c = y_c + l/g*ddy_c;
figure(12); hold on;
plot(t_c, p_c, '-.', 'linewidth', lw, 'color', [0.8500 0.3250 0.0980]);

%% run the simulation with constraints!!

%th = -ddy/g;
%dth = -dddy/g;

%generate the control action
u = griddedInterpolant(t_c, F_c);
data.u = u;

%state is [position, angle] --- [x, xDot, th, thDot]
x0 = [p0; v0; th0; w0];

%run and visualize the simulation
[t_sim, p_sim, v_sim, th_sim, w_sim] = ...
                    Simulation(x0, linspace(t0, tf, nt), data);



%%% Verify we satisfy all of the boundary/interior constraints
figure(10); clf; hold on;
plot(t, y_c, 'linewidth', lw);
plot(t, dy_c, 'linewidth', lw);
plot(t, ddy_c, 'linewidth', lw);
plot(t, dddy_c, 'linewidth', lw);
legend('y', 'y^{(1)}', 'y^{(2)}', 'y^{(3)}');

%check the tangency conditions
p_test = ya(end) + l/g*ddya(end);
v_test = dya(end) + l/g*dddya(end);
fprintf("N1: %g, wanted %g\n", p_test, p_max);
fprintf("N2: %g, wanted %g\n", v_test, 0);


p = y_c + l/g*ddy_c;
theta = -ddy_c/g;

if VIZ
    figure(2)
    Visualize_SmallAngle(y_c, ddy_c, t_c, data, 2);
    figure(3);
    Visualize(p_sim, th_sim, t_sim, data, 3);
end

%%%% make figure 12 pretty
figure(12); hold on;
plot(t_sim, p_sim, ':', 'linewidth', lw, 'color', [0.4660 0.6740 0.1880])
legend('Bound', 'Unconstrained', 'Constrained Small Angle', 'Constrained');

xlabel('Time (s)');
ylabel('p (m)');

grid on; box on;
set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times')


%%%% Make figure 11 pretty
figure(11); hold on;
plot(t_sim, th_sim, ':', 'linewidth', lw, 'color', [0.4660 0.6740 0.1880])
legend('Unconstrained', 'Constrained Small Angle', 'Constrained');

axis([0, tf, -.26, .26])



xlabel('Time (s)');
ylabel('\theta (rad)');

grid on; box on;
set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Times')





%enforce x1 <= x2 -> x1 - x2 <= 0
function [c,ceq] = const(x)
    c = x(1) - x(2);
    ceq = [];
end