function result = open_loop(obj)
%OPEN_LOOP

nx = obj.nx;
N  = obj.ctrl.N;
result = struct;

% number of points in the integration 
% Nfull = 3 * obj.ctrl.Nsub * N;

% initial condition
x0 = obj.output.x(1:nx);

% time span
% if (obj.initial_time_free && obj.final_time_free)
%     t0 = obj.output.p(2);
%     tf = obj.output.p(1);
% elseif (obj.final_time_free)
%     t0 = obj.bnds.init.t_min;
%     tf = obj.output.p(1);
% elseif (obj.initial_time_free)
%     t0 = obj.output.p(1);
%     tf = obj.bnds.trgt.t_min;
% end
tspan = get_time(obj);

% get control matrix (reshape to nu by N) and parameters p
u   = reshape( obj.output.u, obj.nu, N );
% ut  = linspace(t0,tf,N);
p   = obj.output.p;

% propagate dynamics
[~,x] = ode45(@(t,x)deriv(t,x,u,p),tspan,x0);
x = x.';
% compute final state error
xf  = obj.output.x(nx*(N-1)+(1:nx));
err = norm(xf-x(:,end),2);

% output of the open loop propagation
result.x        = x;
result.t        = tspan;
result.err      = err;
% result.Nfull    = Nfull;

    function dx = deriv(t,x,u,p)
        kp = obj.auxdata.kp;
        kd = obj.auxdata.kd;
        k = find(t >= tspan, 1, "last");
        if k == 1
            p_idx = 1;
        else
            p_idx = sum(k >= kp) + 1;
        end

        % interpolate control
        if p_idx < obj.np && k == kp(p_idx) - 1 && kd(p_idx) == 0 
            u_t = interp1(tspan,u',t,'previous')';
        else
            u_t = interp1(tspan,u',t,'linear')';
        end
        
        % query dynamics
        dx = obj.dynamics(obj,t,x,u_t,p);
        
    end

end

