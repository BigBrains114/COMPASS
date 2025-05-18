function [A,B,C] = multiphase_linearize(input,t,x,u,p)
%LINEARIZE    
%
% Returns linearized matrices at point (t,x,u,p) for the planar landing
% problem.
%
% Syntax is [A,B,C] = obj.linearize(obj,t,x,u,p) where obj is an scvx
% object and outputs are
%   - A : Jacobian dfdx
%   - B : Jacobian dfdu
%   - C : Jacobian dfdp
%
% T. Reynolds -- RAIN Lab

G = u(1);

if G < 1e-3
    A = A0(x,u);
    B = B0(x,u);
else
    A = A1(x,u);
    B = B1(x,u);
end

C = input.dynamics(input,t,x,u,p);

end

