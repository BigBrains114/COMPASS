function discretize(obj)

    nx  = obj.nx;
    nu  = obj.nu;
    np  = obj.np;
    N   = obj.ctrl.N;
    kp  = obj.auxdata.k_phase;
    kd  = obj.auxdata.k_disc;
    tau = obj.auxdata.tau;

    % Determine if first iteration to use initial guess
    first_iter = isempty(obj.output) || ~isfield(obj.output,'x') || isempty(obj.output.x);

    % Get reference parameters (phase durations)
    if (first_iter)
        pref = obj.iguess.p; % Size np x 1
    else
        pref = obj.output.p;
    end

    % --- Initialize Output Matrices (Sparse Assembly Format) ---
    % Use zeros instead of spalloc, similar to original discretize.m
    EH = zeros(nx*N, nx*N); % Will hold Ad blocks
    BE = zeros(nx*N, nu*N); % Will hold Bm, Bp blocks
    ES = zeros(nx*N, np);   % Will hold S blocks (mapped to parameters)
    AR = zeros(nx*N, 1);   % Will hold R (residual/offset d_k) blocks

    % Set identity for first block of EH (x_1 = x_1)
    % Use scvx.set_block [cite: 5]
    EH = scvx.set_block(EH, eye(nx), 1, 1);

    obj.defects = zeros(N-1, 1); % Initialize defects vector

    % --- Define Augmented State Size and Indices ---
    dim_P = nx + nx*nx + 2*nx*nu + nx; % x, Phi, Psi_Bm, Psi_Bp, Psi_S
    idx_x = 1:nx;
    idx_PA = nx + (1:nx*nx);
    idx_PBm = idx_PA(end) + (1:nx*nu);
    idx_PBp = idx_PBm(end) + (1:nx*nu);
    idx_PS = idx_PBp(end) + (1:nx); % Sensitivity w.r.t. interval duration sk

    % --- Set Initial Conditions for Augmented State's Sensitivity Part ---
    P0 = zeros(dim_P, 1); % Initialize sensitivity part
    P0(idx_PA) = reshape(eye(nx), [], 1); % Phi(0) = Identity

    % --- Loop Through Time Intervals ---
    p_idx = 1; % Current phase index (1 to np)

    for k = 1:N-1
        tspan = linspace(tau(k),tau(k+1),obj.ctrl.Nsub);
        if (first_iter)
            xk = obj.iguess.x(:,k);
            xk1 = obj.iguess.x(:,k+1);
            u  = obj.iguess.u(:,k:k+1);
        else
            xk = scvx.get_block(obj.output.x,k,1,N,1);
            xk1 = scvx.get_block(obj.output.x,k+1,1,N,1);
            u = [ scvx.get_block(obj.output.u,k,1,N,1), ...
                  scvx.get_block(obj.output.u,k+1,1,N,1) ];
        end
        
        % set inital condition for propagated state
        P0(idx_x) = xk;

        % Determine the phase index for the current interval k
        if p_idx < np && k > kp(p_idx) - 1
             p_idx = p_idx + 1;
        end
        
        % get time-dilation factor for current phase
        sk = pref(p_idx);

        % integrate 
        P = scvx.rk4(@(t,P)deriv(t,P,u,sk),tspan,P0);
        P = P(end, :); % Get the final augmented state vector

        xpk1 = reshape(P(idx_x),nx,1);
        Ak   = reshape(P(idx_PA), nx, nx);
        Bmk  = reshape(P(idx_PBm), nx, nu);
        if p_idx < np && k == kp(p_idx) - 1 && kd(p_idx) == 0 
            Bpk   = zeros(nx, nu);
        else
            Bpk   = reshape(P(idx_PBp), nx, nu);
        end
        Sk_int    = reshape(P(idx_PS),nx,1);

        % --- Calculate Residual Term d_k (Eq. 11 format) ---
        % d_k = x_prop_{k+1} - (A_k*x_k_ref + Bm_k*u_k_ref + Bp_k*u_{k+1}_ref + Sk_int*sk_ref)
        dk = xpk1 - (Ak * xk + Bmk * u(:,1) + Bpk * u(:,2) + Sk_int * sk);

        % --- Assemble Sparse Matrices using scvx.set_block [cite: 5] ---
        % Note: Row index k+1 corresponds to constraint linking k and k+1
        EH = scvx.set_block(EH, Ak,      k+1, k);     % Coefficient of x_k
        % EH already has Identity on diagonal from initialization or previous step

        BE = scvx.set_block(BE, Bmk,     k+1, k);     % Coeff of u_k
        BE = scvx.set_block(BE, Bpk,     k+1, k+1);   % Coeff of u_{k+1}

        % ES: Contains Sk mapped to the correct parameter column
        ES = scvx.set_block(ES, Sk_int,  k+1, p_idx); % Coeff of p_phase_idx

        % AR: Contains the residual/offset term dk
        AR = scvx.set_block(AR, dk,      k+1, 1);

        % --- Optional: Calculate Defect (for monitoring) ---
        defect = norm(xk1 - xpk1, 2);
        obj.defects(k) = defect;

    end

    % --- Store Results in obj.output ---
    % Ad includes the negative of A_k blocks AND the Identity blocks
    % Need to subtract Identity to match original Ad definition if needed,
    % but EH as constructed here represents: I*x_{k+1} - Ak*x_k = ...
    obj.output.Ad   = EH; % Directly use EH (includes I blocks)
    obj.output.Bd   = BE; % Directly use BE (includes Bm, Bp blocks)
    obj.output.Sd   = ES;
    obj.output.Rd   = AR; % Corresponds to the offset d_k

    % --- Assess Feasibility (Optional, based on defect) ---
    feas = all(obj.defects <= obj.ctrl.feas_tol);
    obj.output.feas = feas;


    function dP = deriv(tau,P,u,p)
        x_tau = reshape(P(idx_x), nx, 1);
        PA    = reshape(P(idx_PA), nx, nx);
        PBm   = reshape(P(idx_PBm), nx, nu);
        PS    = reshape(P(idx_PS), nx, 1);

        if p_idx < np && k == kp(p_idx) - 1 && kd(p_idx) == 0
            u_tau = interp1([tspan(1), tspan(end)],u',tau,'previous')';
            
            f = obj.dynamics(obj, tau, x_tau, u_tau, p);
            [A,B,C] = obj.linearize(obj,tau,x_tau,u_tau,p);
            
            dBm  = (p * A * PBm) + (p * B);
            dBp  = zeros(nx, nu);
        else
            u_tau = interp1([tspan(1), tspan(end)],u',tau,'linear')';
            PBp   = reshape(P(idx_PBp), nx, nu);
            
            lm    = (tspan(end)-tau)/(tspan(end)-tspan(1));
            lp    = (tau-tspan(1))/(tspan(end)-tspan(1));
            
            f = obj.dynamics(obj, tau, x_tau, u_tau, p);
            [A,B,C] = obj.linearize(obj,tau,x_tau,u_tau,p);

            dBm  = (p * A * PBm) + (p * B * lm);
            dBp  = (p * A * PBp) + (p * B * lp);
        end

        dx = p * f;
        dPA = (p * A) * PA;
        dS   = (p * A * PS) + C;

        dP = [ dx; dPA(:); dBm(:); dBp(:); dS(:) ];
    end 

end 