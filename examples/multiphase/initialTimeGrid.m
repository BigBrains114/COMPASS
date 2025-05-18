function t = initialTimeGrid(obj)
    N = obj.ctrl.N;
    np = obj.np;
    kp = obj.auxdata.kp;
    x0 = obj.auxdata.x0;
    
    % if k_phase(1) ~= 1
    %     k_phase = [1 k_phase];
    % end

    % solve analytic equation for time of flight
    tf = roots([9.81^2/2, 0, -2*dot(x0(4:5),x0(4:5)), -12*dot(x0(4:5),x0(2:3)), -18*dot(x0(2:3),x0(2:3))]);
    tf = tf(~imag(tf) & real(tf)>0); % choose the correct solution

    % Calculate the duration of each phase (divide tf evenly across phases)
    phase_durations = tf / np * ones(np, 1);
    num_points = diff([1; kp(:); N]);
    s = phase_durations ./ num_points;
    
    % Initialize time vector
    t = zeros(1, N);
    kp = [1 kp];
    % sk = s(1);
    for k = 1:N-1
        % if any(k >= k_phase)
        sk = s(find(k >= kp,1,"last"));
        % end
        t(k+1) = t(k) + sk;
    end

    % last_time = t(1);
    % start_idx = 1;
    % Distribute points in each phase
    % for id = 1:np
    %     if id < np
    %         % end_idx = k_phase(id+1) - 1;
    %         % num_delta = k_phase(id) - start_idx;
    %         t(k_phase(id)) = last_time + phase_durations(id);
    % 
    %         if num_delta > 1
    %             t_grid = linspace(last_time,t(k_phase(id)),num_delta+1);
    %             t(start_idx+1:k_phase(id)-1) = t_grid(2:end-1);
    %         end
    % 
    %         last_time = t(k_phase(id));
    %     else
    %         % end_idx = N; % For the last phase
    %         % num_delta = N - k_phase(id) + 1;
    %         t(end) = last_time + phase_durations(id);
    %     end
    % 
        % start_idx = k_phase(id);
        % Start and end times of the current phase
        % if id == 1
        %     t_start = 0;
        % else
        %     t_start = sum(phase_durations(1:id-1));
        % end
        % t_end = t_start + phase_durations(id);
        % 
        % if num_delta > 1
        %     % Create evenly spaced time points for the current phase
        %     t((start_idx+1):end_idx) = linspace(t_start, t_end, num_delta);
        % else
        %     t(start_idx) = t_start;
        % end

        
    % end
end