function controls = initialControlGuess(obj, x0)
    N = obj.ctrl.N;
    nu = obj.nu;
    kp = obj.auxdata.kp;
    % g0 = 9.81;
    % Tmax = obj.bnds.u_max(1);
    % Tmin = obj.bnds.u_min(1);

    Tmax = 2200e3;
    Tmin = 880e3;

    controls = zeros(nu, N);
    % for id = 1:N
    %     controls(1,id) = x0(1,id)*g0;
    %     if controls(1,id) > Tmax
    %         controls(1,id) = Tmax;
    %     elseif controls(1,id) < Tmin
    %         controls(1,id) = Tmin;
    %     end
    % end

    controls(:,1)               = [0;0];
    controls(:,kp(1):kp(2)-1)   = repmat([3*Tmax;  0],1,kp(2)-kp(1));
    controls(:,kp(2):kp(3)-1)   = repmat([  Tmax;  0],1,kp(3)-kp(2));
    controls(:,kp(3):N      )   = repmat([  Tmin;  0],1,  N-kp(3)+1);

    % controls(controls > Tmax) = Tmax;
    % controls(controls < Tmin) = Tmin;
end