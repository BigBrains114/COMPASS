function plot_trajectory_data(obj, result_open_loop, iter)
% plot_trajectory_data Generates plots for states and controls from an scvx object.
%
% Inputs:
%   obj: The scvx object (e.g., 'O' from your main.m script) containing
%        the optimization output in obj.output, and problem definitions
%        in obj.nx, obj.nu, obj.ctrl, obj.bnds.
%   result_open_loop: (Optional) The result structure from an open-loop
%                     simulation (e.g., from obj.open_loop()). If provided,
%                     the propagated trajectory will be plotted alongside
%                     the discrete solution. If empty or not provided,
%                     only the discrete solution is plotted.

% --- Input Handling ---
if nargin < 2
    result_open_loop = []; % Set to empty if not provided
end
plot_open_loop = ~isempty(result_open_loop) && isstruct(result_open_loop) && ...
                 isfield(result_open_loop, 't') && isfield(result_open_loop, 'x');

% --- Setup Plotting Style ---
set(0,'defaulttextinterpreter','latex','defaultAxesFontSize',16,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex',...
    'defaultLineLineWidth',1.5,...
    'defaultLineMarkerSize',6); % Increased marker size slightly

% --- Extract Data from scvx Object ---
% if ~isfield(obj, 'output') || ~isfield(obj.output, 'x') || ~isfield(obj.output, 'u')
%     error('plot_trajectory_data:MissingOutput', ...
%           'The scvx object does not contain the required output data (obj.output.x, obj.output.u).');
% end

X_discrete = reshape(obj.output.x, obj.nx, obj.ctrl.N);
U_discrete = reshape(obj.output.u, obj.nu, obj.ctrl.N);
T_discrete = get_time(obj); % Assuming get_time is a method of the scvx object or accessible

% State variables based on your definition:
% x = [mass; altitude; horizontal_pos; altitude_rate; horizontal_pos_rate; attitude; attitude_rate]
m_discrete  = X_discrete(1,:);
pos_discrete = X_discrete(2:3,:); % rows: altitude, horizontal_pos
vel_discrete = X_discrete(4:5,:); % rows: altitude_rate, horizontal_pos_rate
att_discrete = X_discrete(6,:);
ang_rate_discrete = X_discrete(7,:);

% Time limits for plots
Tlim_discrete = [0 T_discrete(end)];
if T_discrete(end) == 0 && obj.ctrl.N > 1 % Handle case of zero total time if all p are zero
    Tlim_discrete = [0 1]; % Default if time is zero
end

% --- Figure 1: Mass, Position, Velocity ---
f1 = figure(1);
clf;
f1.Name = 'Mass, Position, Velocity vs. Time';
% Mass
subplot(3,1,1); hold on; grid on; box on;
plot(T_discrete, m_discrete, 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimized Mass');
if plot_open_loop
    plot(result_open_loop.t, result_open_loop.x(1,:), 'b-', 'DisplayName', 'Open-Loop Mass');
end
set(gca, 'Xlim', Tlim_discrete);
ylabel('Mass ($kg$)');
title('Mass vs. Time');
legend('show', 'Location', 'best');

% Position (Altitude and Horizontal Position)
subplot(3,1,2); hold on; grid on; box on;
plot(T_discrete, pos_discrete(1,:), 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Opt. Altitude');
plot(T_discrete, pos_discrete(2,:), 'ro-', 'MarkerFaceColor', 'r', 'DisplayName', 'Opt. Horiz. Pos.');
if plot_open_loop
    plot(result_open_loop.t, result_open_loop.x(2,:), 'b-', 'DisplayName', 'OL Altitude');
    plot(result_open_loop.t, result_open_loop.x(3,:), 'm-', 'DisplayName', 'OL Horiz. Pos.');
end
set(gca, 'Xlim', Tlim_discrete);
ylabel('Position ($m$)');
title('Position vs. Time');
legend('show', 'Location', 'best');

% Velocity (Altitude Rate and Horizontal Position Rate)
subplot(3,1,3); hold on; grid on; box on;
plot(T_discrete, vel_discrete(1,:), 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Opt. Altitude Rate');
plot(T_discrete, vel_discrete(2,:), 'ro-', 'MarkerFaceColor', 'r', 'DisplayName', 'Opt. Horiz. Rate');
if plot_open_loop
    plot(result_open_loop.t, result_open_loop.x(4,:), 'b-', 'DisplayName', 'OL Altitude Rate');
    plot(result_open_loop.t, result_open_loop.x(5,:), 'm-', 'DisplayName', 'OL Horiz. Rate');
end
set(gca, 'Xlim', Tlim_discrete);
ylabel('Velocity ($m/s$)');
xlabel('Time ($s$)');
title('Velocity vs. Time');
legend('show', 'Location', 'best');

% --- Figure 2: Attitude and Angular Rate ---
f2 = figure(2);
clf;
f2.Name = 'Attitude, Angular Rate vs. Time';
% Attitude Angle
subplot(2,1,1); hold on; grid on; box on;
plot(T_discrete, rad2deg(att_discrete), 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimized Attitude');
if plot_open_loop
    plot(result_open_loop.t, rad2deg(result_open_loop.x(6,:)), 'b-', 'DisplayName', 'Open-Loop Attitude');
end
set(gca, 'Xlim', Tlim_discrete);
ylabel('Attitude Angle ($deg$)');
title('Attitude Angle vs. Time');
legend('show', 'Location', 'best');

% Angular Rate
subplot(2,1,2); hold on; grid on; box on;
plot(T_discrete, rad2deg(ang_rate_discrete), 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimized Ang. Rate');
if plot_open_loop
    plot(result_open_loop.t, rad2deg(result_open_loop.x(7,:)), 'b-', 'DisplayName', 'Open-Loop Ang. Rate');
end
set(gca, 'Xlim', Tlim_discrete);
ylabel('Angular Rate ($deg/s$)');
xlabel('Time ($s$)');
title('Angular Rate vs. Time');
legend('show', 'Location', 'best');

% --- Figure 3: Controls ---
f3 = figure(3);
clf;
f3.Name = 'Controls vs. Time';
% Thrust
subplot(2,1,1); hold on; grid on; box on;
plot(T_discrete, U_discrete(1,:), 'bo-', 'MarkerFaceColor', 'b', 'DisplayName', 'Optimized Thrust');
if isfield(obj, 'bnds') && isfield(obj.bnds, 'u_min') && isfield(obj.bnds, 'u_max')
    plot(Tlim_discrete, [obj.bnds.u_min(1) obj.bnds.u_min(1)], 'r--', 'DisplayName', 'Min Thrust Bound');
    plot(Tlim_discrete, [obj.bnds.u_max(1) obj.bnds.u_max(1)], 'r--', 'DisplayName', 'Max Thrust Bound');
    % For multi-phase thrust limits, this simple bound plot might not capture all details.
    % You might need to plot piecewise bounds if Tmin/Tmax change per phase.
end
set(gca, 'Xlim', Tlim_discrete);
if isfield(obj, 'bnds') && isfield(obj.bnds, 'u_max')
    ylim_padding = (obj.bnds.u_max(1) - obj.bnds.u_min(1)) * 0.1;
    if ylim_padding == 0; ylim_padding = 1; end % Avoid zero padding
    set(gca, 'Ylim', [obj.bnds.u_min(1) - ylim_padding, obj.bnds.u_max(1) + ylim_padding]);
end
ylabel('Thrust ($N$)');
title('Thrust vs. Time');
legend('show', 'Location', 'best');

% Gimbal Angle
subplot(2,1,2); hold on; grid on; box on;
plot(T_discrete, rad2deg(U_discrete(2,:)), 'bo-', 'MarkerFaceColor', 'b', 'DisplayName', 'Optimized Gimbal Angle');
if isfield(obj, 'bnds') && isfield(obj.bnds, 'u_min') && isfield(obj.bnds, 'u_max')
    plot(Tlim_discrete, rad2deg([obj.bnds.u_min(2) obj.bnds.u_min(2)]), 'r--', 'DisplayName', 'Min Gimbal Bound');
    plot(Tlim_discrete, rad2deg([obj.bnds.u_max(2) obj.bnds.u_max(2)]), 'r--', 'DisplayName', 'Max Gimbal Bound');
end
set(gca, 'Xlim', Tlim_discrete);
if isfield(obj, 'bnds') && isfield(obj.bnds, 'u_min') && isfield(obj.bnds, 'u_max')
     ylim_padding_deg = rad2deg(obj.bnds.u_max(2) - obj.bnds.u_min(2)) * 0.1;
     if ylim_padding_deg == 0; ylim_padding_deg = 5; end % Avoid zero padding
    set(gca, 'Ylim', rad2deg([obj.bnds.u_min(2) - deg2rad(ylim_padding_deg/2), obj.bnds.u_max(2) + deg2rad(ylim_padding_deg/2)]));
end
ylabel('Gimbal Angle ($deg$)');
xlabel('Time ($s$)');
title('Gimbal Angle vs. Time');
legend('show', 'Location', 'best');

% --- Figure 4: Trajectory ---
f4 = figure(4);
clf; hold on; grid on; box on;
f4.Name = 'Landing Trajectory';
plot(pos_discrete(2,:), pos_discrete(1,:), 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Optimized Trajectory');
if plot_open_loop
    plot(result_open_loop.x(3,:), result_open_loop.x(2,:), 'b-', 'DisplayName', 'Open-Loop Trajectory');
end
% Plot target landing spot
if isfield(obj, 'auxdata') && isfield(obj.auxdata, 'xf') && length(obj.auxdata.xf) >= 2
    % Assuming xf(1) is target altitude, xf(2) is target horizontal position
    % based on xf_idx = 2:7. So auxdata.xf(1) is target for state 2 (altitude),
    % auxdata.xf(2) is target for state 3 (horizontal position).
    plot(obj.auxdata.xf(2), obj.auxdata.xf(1), 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Target');
end
% Plot initial position
plot(pos_discrete(2,1), pos_discrete(1,1), 'go', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');

ylabel('Altitude ($m$)');
xlabel('Horizontal Position ($m$)');
title('Landing Trajectory (Altitude vs. Horizontal Position)');
axis equal; % Important for trajectory plots
legend('show', 'Location', 'best');
set(gca,'YDir','normal'); % Ensure altitude increases upwards

% --- Figure 5: Virtual States (if they exist and are different) ---
if isfield(obj.output, 'vs') && ~isempty(obj.output.vs)
    VS_discrete = reshape(obj.output.vs, obj.nx, obj.ctrl.N);
    
    % Check if virtual states are significantly different from actual states
    if norm(X_discrete - VS_discrete, 'fro') > 1e-3 % Tolerance for difference
        
        f5 = figure(5);
        clf;
        f5.Name = 'Virtual States vs. Actual States';
        
        state_names = {'Mass', 'Altitude', 'Horizontal Pos', 'Altitude Rate', 'Horiz. Pos. Rate', 'Attitude (deg)', 'Ang. Rate (deg/s)'};
        
        for i = 1:obj.nx
            subplot(ceil(obj.nx/2), 2, i); hold on; grid on; box on;
            
            current_X = X_discrete(i,:);
            current_VS = VS_discrete(i,:);
            
            if i == 6 || i == 7 % Convert angles to degrees for plotting
                current_X = rad2deg(current_X);
                current_VS = rad2deg(current_VS);
            end

            plot(T_discrete, current_X, 'ko-', 'MarkerFaceColor', 'k', 'DisplayName', 'Actual State');
            plot(T_discrete, current_VS, 'ro--', 'MarkerFaceColor', 'r', 'DisplayName', 'Virtual State');
            
            set(gca, 'Xlim', Tlim_discrete);
            ylabel(state_names{i});
            title(['State: ', state_names{i}]);
            legend('show', 'Location', 'best');
            if i >= obj.nx - 1 % Add xlabel to bottom plots
                 xlabel('Time ($s$)');
            end
        end
    else
        fprintf('Virtual states are very close to actual states. Skipping virtual state plot.\n');
    end
end


% fprintf('Plotting complete.\n');
set(0,'defaulttextinterpreter','tex'); % Reset to default
end
