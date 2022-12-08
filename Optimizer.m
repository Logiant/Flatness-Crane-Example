function J = Optimizer(x, t0, tf, ICs, FCs, params)
%extract optimization variables t1 and t2
t1 = x(1);
t2 = x(2);

l = params.l;
g = params.g;

ag = params.alpha *g;
w = sqrt(g/l);


%generate the coefficients for this trajectory
c = ConstrainedBVP(t1, t2, t0, tf, ICs, FCs, params);


%calculate how close we are to satisfying continuity in y4
J1 = abs( c(7)/(ag)^4*exp(t1/ag) + c(8)/(ag)^4*exp(-t1/ag) + 24*c(5) + 120*c(6)*t1 - ...
          w^4*( c(9)*cos(w*t1) + c(10)*sin(w*t1)) );

J2 = abs( c(17)/(ag)^4*exp(t2/ag) + c(18)/(ag)^4*exp(-t2/ag) + 24*c(15) + 120*c(16)*t2 - ...
          w^4*( c(9)*cos(w*t2) + c(10)*sin(w*t2)) );

J = J1^2 + J2^2;
      
end

