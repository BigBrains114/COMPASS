function varxu = solve_socp(obj,scaling,iter)
%SOLVE_SOCP
%
% T. Reynolds

nx = obj.nx;
nu = obj.nu;
np = obj.np;
N  = obj.ctrl.N;
kp = obj.auxdata.kp;
xf_idx = obj.auxdata.xf_idx;

first   = 1:nx;
last    = nx*(N-1)+1:nx*N;
% last    = last(xf_idx);
% isptr   = strcmp(obj.ctrl.algo,'ptr');

% reference values
wc      = obj.ctrl.wc;
wvse    = obj.ctrl.wvse;
wtr     = obj.ctrl.wtr;
% wtrp    = obj.ctrl.wtrp;
% if (~isptr)
%     tr = obj.output.tr;
% end
if (iter < 2) 
    xref = obj.iguess.x(:);
    uref = obj.iguess.u(:);
    pref = obj.iguess.p;
    vref = obj.iguess.x(:);
else
    xref = obj.output.x;
    uref = obj.output.u;
    pref = obj.output.p;
    vref = obj.output.vs;
end
EH = obj.output.Ad;
BE = obj.output.Bd;
ES = obj.output.Sd;
AR = obj.output.Rd;

% scaling matrices
Sx = scaling.Sx;
cx = scaling.cx;
SX = scaling.SX;
cX = scaling.cX;
Su = scaling.Su;
cu = scaling.cu;
SU = scaling.SU;
cU = scaling.cU;
Sp = scaling.Sp;
cp = scaling.cp;
iSx = scaling.iSx;
iSu = scaling.iSu;
iSp = scaling.iSp;

% constraint things
% x_min = obj.bnds.x_min;
% x_max = obj.bnds.x_max;
% u_min = obj.bnds.u_min;
% u_max = obj.bnds.u_max;
% nx_constraints_path = numel(obj.bnds.path.x);
% nu_constraints_path = numel(obj.bnds.path.u);

xref_s = zeros(size(xref));
uref_s = zeros(size(uref));
for k = 1:N
   xref_s(nx*(k-1)+(1:nx)) = iSx*(xref(nx*(k-1)+(1:nx))-cx);
   uref_s(nu*(k-1)+(1:nu)) = iSu*(uref(nu*(k-1)+(1:nu))-cu);
end
pref_s = iSp*(pref-cp);


cvx_tic;
cvx_begin quiet
    cvx_solver(obj.ctrl.solver)
    cvx_precision('low')
    
    variable xb(nx*N,1) nonnegative
    variable ub(nu*N,1) nonnegative
    variable pb(np,1) nonnegative   % even in np=0, this is just an empty vector
    variable vs(nx*N,1) nonnegative % changed virtual controls to virtual states
    
    Jtr  = (xb-xref_s).'*(xb-xref_s) + (ub-uref_s).'*(ub-uref_s) + (pb-pref_s).'*(pb-pref_s);
    Jvse = (xb-vs).'*(xb-vs); 
    cost = wc * obj.cost(xb,ub,pb,N) + wtr * Jtr + wvse * Jvse;
    minimize( cost )
    
    subject to
    
    % initial conditions
    % vs(first) == iSx*(obj.auxdata.x0-cx);
    (Sx*vs(first)+cx) == obj.auxdata.x0;
    % final conditions
    % vs_f = (Sx*vs(last)+cx);
    % xf_s = iSx*([0;obj.auxdata.xf]-cx);
    vs_f = vs(last);
    % vs_f(xf_idx) == xf_s(xf_idx);
    (Sx(xf_idx,xf_idx)*vs_f(xf_idx)+cx(xf_idx)) == obj.auxdata.xf;
    % dynamics
    (SX*xb+cX) == EH*(SX*xb+cX) + BE*(SU*ub+cU) + ES*(Sp*pb+cp) + AR;
    
    % upper limit on states and virtual states
    xb <= 1;
    ub <= 1;
    pb <= 1;
    % vs <= 1;
    
    p_idx = 1;
    % time loop
    for k = 1:N
       
        if p_idx < obj.np && k > obj.auxdata.kp(p_idx) - 1
             p_idx = p_idx + 1;
        end
        
        % xk   = xb(nx*(k-1)+1:nx*k);
        uk   = Su*ub(nu*(k-1)+1:nu*k)+cu;
        if k > 1
            ubkm1 = uref(nu*(k-2)+1:nu*(k-1));
        end
        
        vsk   = Sx * vs(nx*(k-1)+1:nx*k) + cx;
        vsbk  = vref(nx*(k-1)+1:nx*k);

        if any(k == kp)
            sbkm1  = pref(p_idx-1);
        else
            sbkm1  = pref(p_idx);
        end
        

        if k >= 1 && k < kp(1)        % unpowered
            
            vsk(2) >= obj.auxdata.h_trig;  % minimum altitude before TD
            uk(1) == 0;           % zero inital thrust

        elseif k >= kp(1) && k < kp(2)    % high thrust
            
            vsk(2) >= obj.auxdata.h_trig;  % minimum altitude before TD
            
            % control bounds (eq 24 and 25) 
            uk >= [3*obj.auxdata.Tmin; max(obj.bnds.u_min(2), -obj.auxdata.dmax*sbkm1+ubkm1(2));];
            uk <= [3*obj.auxdata.Tmax; min(obj.bnds.u_max(2),  obj.auxdata.dmax*sbkm1+ubkm1(2));];

        elseif k >= kp(2) && k < kp(3)    % low thrust

            vsk(2) >= obj.auxdata.h_trig; % minimum altitude before TD

            if k == kp(2)
                uk(2) == ubkm1(2);
                uk(1) >= obj.auxdata.Tmin;
                uk(1) <= obj.auxdata.Tmax;
            else
                % control bounds (eq 27 and 28) 
            uk >= [obj.auxdata.Tmin; max(obj.bnds.u_min(2), -obj.auxdata.dmax*sbkm1+ubkm1(2));];
            uk <= [obj.auxdata.Tmax; min(obj.bnds.u_max(2),  obj.auxdata.dmax*sbkm1+ubkm1(2));];
            end
        elseif k >= kp(3)                 % terminal
            % control bounds (eq 29 and 30f) 
            uk >= [obj.auxdata.Tmin; max(obj.bnds.u_min(2), -obj.auxdata.dmaxTD*sbkm1+ubkm1(2));];
            uk <= [obj.auxdata.Tmax; min(obj.bnds.u_max(2),  obj.auxdata.dmaxTD*sbkm1+ubkm1(2));];
            
            if k == kp
                vsk(2) == obj.auxdata.h_trig;
                abs(vsk(3)) <= tan(obj.auxdata.gs)*obj.auxdata.h_trig;
            else
                vsk(2) <= obj.auxdata.h_trig;
                abs(vsk(3)) <= tan(obj.auxdata.gs)*vsbk(2);
            end

           norm(vsk(4:5)) <= obj.auxdata.vmax;
           abs(vsk(6)) <= obj.auxdata.amax;
           abs(vsk(7)) <= obj.auxdata.wmax;
        else 
            error("Grid point out of range.");
        end

        % linear bounds
        % (Sx*xk+cx) <= x_max;
        % (Sx*xk+cx) >= x_min;
        % (Su*uk+cu) <= u_max;
        % (Su*uk+cu) >= u_min;
        % (Sx*vsk+cx) <= x_max;
        % (Sx*vsk+cx) >= x_min;
        
        % state path constraints
        % for xcnstr = 1:nx_constraints_path
        %     if (obj.bnds.path.x_cvx{xcnstr})
        %         obj.bnds.path.x{xcnstr}( Sx*xbk+cx ) <= 0.0;
        %     else
        %         TS = obj.bnds.path.x_lin{xcnstr};
        %         TS{1}(k) + TS{2}(k,:) * (Sx*xbk+cx) <= 0.0;
        %     end
        % end
                       
        % control path constraints
        % for ucnstr = 1:nu_constraints_path
        %     if (obj.bnds.path.u_cvx{ucnstr})
        %         obj.bnds.path.u{ucnstr}( Su*ubk+cu ) <= 0.0;
        %     else
        %         TS = obj.output.bnds.u{ucnstr};
        %         TS{1} + TS{2} * (Su*ubk+cu) <= 0.0;
        %     end
        % end
        
    end
cvx_end
timing = cvx_toc;

% obj.dbg.out.xb  = xb;
% obj.dbg.out.ub  = ub;
% obj.dbg.out.pb  = pb;
% obj.dbg.out.vs  = vs;


obj.output.x    = SX*xb+cX;
obj.output.u    = SU*ub+cU;
obj.output.p    = Sp*pb+cp;
obj.output.vs   = SX*vs+cX;
obj.output.vse  = Jvse;
obj.output.tr   = Jtr;
obj.output.cost = cvx_optval;


% disp("final state error norm: " + norm(obj.output.x(xf_idx + (N-1)*nx) - obj.auxdata.xf))

% compute max (scaled) change in state and/or control
temp = 0.0;
for k = 1:N
   tempx = norm(xb(nx*(k-1)+(1:nx)) - iSx*(xref(nx*(k-1)+(1:nx))-cx),inf);
   tempu = norm(ub(nu*(k-1)+(1:nu)) - iSu*(uref(nu*(k-1)+(1:nu))-cu),inf);
   temp  = max([temp,tempx,tempu]);    
end
varxu = temp;

% save some data from the solution
if (iter>1)
    obj.output.data.cum_time = obj.output.data.cum_time + sum(timing(4:5));
else
    obj.output.data.cum_time = sum(timing(4:5));
end
switch cvx_status
    case 'Solved'
        obj.output.data.status{iter} = 0;
    case 'Inaccurate/Solved'
        obj.output.data.status{iter} = 1;
    case {'Failed','Unbounded','Inaccurate/Unbounded'}
        obj.output.data.status{iter} = 2;
    case {'Infeasible','Inaccurate/Infeasible'}
        obj.output.data.status{iter} = 3;
    otherwise
        obj.output.data.status{iter} = 4;
end
obj.output.data.slvitr = cvx_slvitr;

end

