function C = UnconstrainedBVP(ICs, FCs, t, params)

alpha = params.alpha;
g = params.g;

t0 = t(1); tf=t(end);

A = [1, t0, t0^2,   t0^3,   t0^4,    t0^5, exp(t0/(alpha*g)),              exp(-t0/(alpha*g)); ...
     0,  1, 2*t0, 3*t0^2, 4*t0^3,  5*t0^4, exp(t0/(alpha*g))/(alpha*g),   -exp(-t0/(alpha*g))/(alpha*g); ...
     0,  0, 2,    6*t0,  12*t0^2, 20*t0^3, exp(t0/(alpha*g))/(alpha*g)^2,  exp(-t0/(alpha*g))/(alpha*g)^2; ...
     0,  0, 0,    6,     24*t0,   60*t0^2, exp(t0/(alpha*g))/(alpha*g)^3, -exp(-t0/(alpha*g))/(alpha*g)^3; ...
     1, tf, tf^2,   tf^3,   tf^4,    tf^5, exp(tf/(alpha*g)),              exp(-tf/(alpha*g)); ...
     0,  1, 2*tf, 3*tf^2, 4*tf^3,  5*tf^4, exp(tf/(alpha*g))/(alpha*g),   -exp(-tf/(alpha*g))/(alpha*g); ...
     0,  0, 2,    6*tf,  12*tf^2, 20*tf^3, exp(tf/(alpha*g))/(alpha*g)^2,  exp(-tf/(alpha*g))/(alpha*g)^2; ...
     0,  0, 0,    6,     24*tf,   60*tf^2, exp(tf/(alpha*g))/(alpha*g)^3, -exp(-tf/(alpha*g))/(alpha*g)^3];

b = [ICs; FCs];
C = A\b;

end
