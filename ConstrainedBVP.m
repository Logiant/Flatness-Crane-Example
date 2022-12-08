function C = ConstrainedBVP(t1, t2, t0, tf, ICs, FCs, params)

l = params.l;
g = params.g;
p_max = params.p_max;


ag = params.alpha *g;
w = sqrt(g/l);
%constraints C are c0 c1 c2 c3 c4 c5 c6 c7 | e1 e2 | f0 f1 f2 f3 f4 f5 f6 f7
% A will be an 18x18 matrix with everything except continuity in y^(4)

%initial conditions (b = ICs here)
A0 = [1, t0, t0^2,   t0^3,   t0^4,    t0^5, exp(t0/ag),              exp(-t0/ag); ...
      0,  1, 2*t0, 3*t0^2, 4*t0^3,  5*t0^4, exp(t0/ag)/(ag),   -exp(-t0/ag)/(ag); ...
      0,  0, 2,    6*t0,  12*t0^2, 20*t0^3, exp(t0/ag)/(ag)^2,  exp(-t0/ag)/(ag)^2; ...
      0,  0, 0,    6,     24*t0,   60*t0^2, exp(t0/ag)/(ag)^3, -exp(-t0/ag)/(ag)^3];

%%% First Junction
%tangency condition 1
A1T1 = [1, t1, t1^2+2*l/g, t1^3+6*t1*l/g, t1^4+12*t1^2*l/g, t1^5+20*l/g*t1^3, exp(t1/ag)*(1 + l/g/(ag)^2), exp(-t1/ag)*(1 + l/g/(ag)^2)];
%tangency condition 2
A1T2 = [0, 1, 2*t1, 3*t1^2+6*l/g, 4*t1^3+24*l/g*t1, 5*t1^4+60*l/g*t1^2, exp(t1/ag)/ag + l/g*exp(t1/ag)/(ag)^3, -exp(-t1/ag)/ag - l/g*exp(-t1/ag)/(ag)^3];
%continuity in y
A1C0 = [1, t1, t1^2,t1^3,t1^4,t1^5,exp(t1/ag),exp(-t1/ag),-cos(w*t1),-sin(w*t1)];
%continuity in yDot
A1C1 = [0, 1, 2*t1,3*t1^2,4*t1^3,5*t1^4,exp(t1/ag)/ag,-exp(-t1/ag)/ag,w*sin(w*t1),-w*cos(w*t1)];
%continuity in yDDot
A1C2 = [0,0,2,6*t1,12*t1^2,20*t1^3,exp(t1/ag)/ag^2,exp(-t1/ag)/ag^2,w^2*cos(w*t1),w^2*sin(w*t1)];
%continuity in yDDDot
A1C3 = [0,0,0,6,24*t1,60*t1^2,exp(t1/ag)/ag^3,-exp(-t1/ag)/ag^3,-w^3*sin(w*t1),w^3*cos(w*t1)];

%stack the A matrices and construct B
A1 = [A1T1, zeros(1, 2);
      A1T2, zeros(1, 2);
      A1C0;
      A1C1;
      A1C2;
      A1C3];
B1 = [p_max; 0; p_max; 0; 0; 0];


%%% Second Junction
%continuity in y
A2C0 = [cos(w*t2),sin(w*t2),-1,-t2,-t2^2,-t2^3,-t2^4,-t2^5,-exp(t2/ag),-exp(-t2/ag)];
%continuity in yDot
A2C1 = [-w*sin(w*t2),w*cos(w*t2),0,-1,-2*t2,-3*t2^2,-4*t2^3,-5*t2^4,-exp(t2/ag)/ag,exp(-t2/ag)/ag];
%continuity in yDDot
A2C2 = [-w^2*cos(w*t2),-w^2*sin(w*t2),0,0,-2,-6*t2,-12*t2^2,-20*t2^3,-exp(t2/ag)/ag^2,-exp(-t2/ag)/ag^2];
%continuity in yDDDot
A2C3 = [w^3*sin(w*t2),-w^3*cos(w*t2),0,0,0,-6,-24*t2,-60*t2^2,-exp(t2/ag)/ag^3,exp(-t2/ag)/ag^3];
 
%stack the A matrices and construct B
A2 = [A2C0; A2C1; A2C2; A2C3];
B2 = [-p_max; 0; 0; 0];
  
%final conditions (b = FCs here)
AF = [1, tf, tf^2,   tf^3,   tf^4,    tf^5, exp(tf/(ag)),              exp(-tf/(ag)); ...
      0,  1, 2*tf, 3*tf^2, 4*tf^3,  5*tf^4, exp(tf/(ag))/(ag),   -exp(-tf/(ag))/(ag); ...
      0,  0, 2,    6*tf,  12*tf^2, 20*tf^3, exp(tf/(ag))/(ag)^2,  exp(-tf/(ag))/(ag)^2; ...
      0,  0, 0,    6,     24*tf,   60*tf^2, exp(tf/(ag))/(ag)^3, -exp(-tf/(ag))/(ag)^3];
  

%%%% CONSTRUCT THE LINEAR MATRIX

B = [ICs; B1; B2; FCs];

A = [A0,           zeros(4,10);
     A1,           zeros(6, 8);
     zeros(4, 8),  A2;
     zeros(4, 10), AF];
    
C = A\B;

end

