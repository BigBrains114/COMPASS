function scaling = set_scales(obj)
%SET_SCALES
%
% T. Reynolds -- RAIN Lab

x_min = obj.bnds.x_min;
x_max = obj.bnds.x_max;
u_min = obj.bnds.u_min;
u_max = obj.bnds.u_max;
p_max = obj.bnds.p_max;
p_min = obj.bnds.p_min;

nx = obj.nx;
nu = obj.nu;
np = obj.np;
N = obj.ctrl.N;

% choose scaled variable intervals
intrvl_x  = [0;1];
wdth_x    = intrvl_x(2)-intrvl_x(1);
intrvl_u  = [0;1];
wdth_u    = intrvl_u(2)-intrvl_u(1);
intrvl_p  = [0;1];
wdth_p    = intrvl_p(2)-intrvl_p(1);

scaling = struct;

% State terms
Sx  = eye(nx);
iSx = eye(nx);
cx  = zeros(nx,1);
for k = 1:nx
    Sx(k,k)     = ( x_max(k)-x_min(k) ) / wdth_x;
    iSx(k,k)    = 1.0/Sx(k,k);
    cx(k)       = x_min(k) - Sx(k,k) * intrvl_x(1);
end
scaling.Sx  = Sx;
scaling.iSx = iSx;
scaling.cx  = cx;
scaling.SX  = kron(eye(N),Sx);
scaling.iSX = kron(eye(N),iSx);
scaling.cX  = repmat(cx,N,1);

% Control terms
Su  = eye(nu);
iSu = eye(nu);
cu  = zeros(nu,1);
for k = 1:nu
    Su(k,k)     = ( u_max(k)-u_min(k) ) / wdth_u;
    iSu(k,k)    = 1.0/Su(k,k);
    cu(k)       = u_min(k) - Su(k,k) * intrvl_u(1);
end
scaling.Su  = Su;
scaling.iSu = iSu;
scaling.cu  = cu;
scaling.SU  = kron(eye(N),Su);
scaling.iSU = kron(eye(N),iSu);
scaling.cU  = repmat(cu,N,1);

% Parameter terms
Sp  = eye(np);
iSp = eye(np);
cp  = zeros(np,1);
for k = 1:np
    Sp(k,k) = ( p_max(k)-p_min(k) ) / wdth_p;
    iSp(k,k) = 1.0/Sp(k,k);
    cp(k) = p_min(k) - Sp(k,k) * intrvl_p(1);
end
scaling.Sp  = Sp;
scaling.cp  = cp;
scaling.iSp = iSp;


end

