function [t_sim, p_sim, v_sim, th_sim, w_sim] = Simulation(x0,t, data)

%construct the helper ODE
ode = @(t, y) ChenODE(t, y, data);

%do the optimization
[t_sim, x] = ode45(ode, t, x0);

p_sim = x(:,1); v_sim = x(:,2); th_sim = x(:,3); w_sim = x(:,4);
end

