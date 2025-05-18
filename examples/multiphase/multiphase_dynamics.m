function dx = multiphase_dynamics(input,t,x,u,p)
%PLANARLANDING_DYNAMICS     
%
% Continuous dynamic equations for the planar landing problem.
%
% Syntax is dx = obj.dynamics(obj,t,x,u,p) where obj is an scvx object.
%
% T. Reynolds -- RAIN Lab

% gravity
g0      = 9.81;         % Gravity Acceleration at sea level [m/s2]
g       = [-g0; 0];
alpha   = 1/(g0*330);

% geometry 
lr      = 4.5;          % fuselage radius [m]
lh      = 50;           % fuselage height [m]
lcm     = 0.4*lh;       % thrust moment-arm [m]

% aerodynamics
lcp_0   = 0;
lcp_1   = 0.2*lh;
rho     = 1.225;
S_aero  = 545;
c_xhat  = 0.0522;
c_zhat  = 0.4068;


m = x(1);
r = x(2:3);
v = x(4:5);
a = x(6);
w = x(7);

T = u(1);
d = u(2);

% Moment of Inertia
J = m * (lr^2/4+lh^2/12);

% rotation matrix
R_B2I = [cos(a), sin(a); -sin(a), cos(a)];

% thrust force
F_i = T .* [cos(a+d); -sin(a+d)];
F_b = T .* [cos(d  ); -sin(d  )];

% aerodynamic force
A_b = -0.5*rho*S_aero*norm(v)*diag([c_xhat, c_zhat])*R_B2I.'*v;
A_i = R_B2I * A_b;

if T < 1e-3
    lcp = lcp_0;
else
    lcp = lcp_1;
end


% dynamics
dm = -alpha .* T;
dr = v;
dv = (F_i + A_i)./m + g;
da = w;
dw = (F_b(2).*lcm - A_b(2).*lcp)./J;

dx = [ dm; dr; dv; da; dw ];

end
