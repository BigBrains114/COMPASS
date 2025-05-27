% initialize example
run('../scvx_startup');

% auxiliary parameters
auxdata         = struct;
auxdata.kp      = [2 7 12];
auxdata.kd      = [0 0 1];
auxdata.x0      = [1e5; 1e3; 1e2; -90; 0; pi/2; 0];
auxdata.xf      = [0; 0; 0; 0; 0; 0];
auxdata.xf_idx  = 2:7; 
auxdata.h_trig  = 100;
auxdata.Tmin    = 880e3;
auxdata.Tmax    = 2200e3;
auxdata.dmax    = deg2rad(10);
auxdata.dmaxTD  = deg2rad(1);
auxdata.gs      = deg2rad(5);
auxdata.vmax    = 20;
auxdata.amax    = deg2rad(5);
auxdata.wmax    = deg2rad(2.5);

% boundary conditions
bnds = struct;
%            [ m; x; z; vx; vz; a; w; ];   
bnds.x_min = [  85e3;   0; -1e2; -1e2; -3e1; -pi; -pi/2 ];
bnds.x_max = [ 100e3; 1e3;  1e2;  1e1;  3e1;  pi;  pi/2 ];
bnds.u_min = [        0; -deg2rad(10) ];
bnds.u_max = [ 3*2200e3;  deg2rad(10) ];
bnds.p_min = [ .6; .6; .6; .6;]; 
bnds.p_max = [ 10; 2; 2; 2;];
bnds.path.x     = {};
bnds.path.x_cvx = {};
bnds.path.u     = {};
bnds.path.u_cvx = {};

ctrl = struct;
ctrl.N          = 16;
ctrl.Nsub       = 10;
ctrl.iter_max   = 15;
ctrl.wc         = 1e0;
ctrl.wvse       = 1e2;
ctrl.wtr        = 1e1;
% ctrl.wtrp       = 1e-2;
ctrl.cvrg_tol   = 1e-1;
ctrl.feas_tol   = 1e0;
% ctrl.rho0       = 0.0;
% ctrl.rho1       = 0.1;
% ctrl.rho2       = 0.9;
% ctrl.tr_lb      = 0.001;
% ctrl.tr_ub      = 10;
% ctrl.alpha      = 2.0;
% ctrl.beta       = 2.0;
ctrl.algo       = 'ptr';
ctrl.solver     = 'ecos';

% initialize scvx problem with a name and dimensions (name,nx,nu,np)
O = scvx('multiphase',7,2,4);

% attach structs to the created object
O.cost      = @(x,u,p,N)( -x(7*(ctrl.N-1)+1) ); % cost function
O.dynamics  = @multiphase_dynamics;        % main dynamics function
O.linearize = @multiphase_linearize;       % main linearization function
O.auxdata   = auxdata;
O.bnds      = bnds;
O.ctrl      = ctrl;
O.time      = initialTimeGrid(O);

% generate initial guess
% x0      = zeros(O.nx,O.ctrl.N);
% x0(1,:) = linspace(bnds.x_max(1),0.5*(bnds.trgt.x_max(1)+bnds.trgt.x_min(1)),ctrl.N);
% x0(2,:) = linspace(6,0,ctrl.N);
% x0(3,:) = linspace(24,1,ctrl.N);
% x0(4,:) = linspace(-4,0,ctrl.N);
% x0(5,:) = linspace(-2,-0.1,ctrl.N);
% u0      = zeros(O.nu,O.ctrl.N);
% u0(1,:) = [ bnds.u_max(1)*ones(1,4), bnds.u_min(1)*ones(1,2) bnds.u_max(1)*ones(1,ctrl.N-6) ];
% p0      = [5.0; 1.0; 1.0; 1.0];

x0      = initialStateGuess(O);
u0      = initialControlGuess(O, x0);
p0      = initialParamGuess(O);

% call initializing function -- parameters come first, they are the only
% mandatory thing to guess.
O.init(p0,x0,u0);

% solve the problem or grab a solution from a file
exitcode = O.solve();
% O.attach_file('../cpp_hp/data/optimal_xup.txt');
%%
% perform open loop propagation on the last control
result = O.open_loop();

%% Plot the results
set(0,'defaulttextinterpreter','latex','defaultAxesFontSize',16,...
    'defaultAxesTickLabelInterpreter','latex',...
    'defaultLegendInterpreter','latex',...
    'defaultLineLineWidth',1.5,...
    'defaultLineMarkerSize',4)

% computed (discrete) solution
X   = reshape( O.output.x, O.nx, O.ctrl.N );
U   = reshape( O.output.u, O.nu, O.ctrl.N );
T   = get_time(O);
m   = X(1,:);
r   = X(2:3,:);
v   = X(4:5,:);
th  = X(6,:);
w   = X(7,:);
Tlim = [ 0 T(end) ];

figure, clf
subplot(3,1,1), hold on, grid on, box on
plot(T,m,'ko','MarkerFaceColor','k','HandleVisibility','off')
plot(result.t,result.x(1,:))
set(gca,'Xlim',Tlim)
ylabel('Mass')
subplot(3,1,2), hold on, grid on, box on
plot(T,r,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(2:3,:))
set(gca,'Xlim',Tlim)
ylabel('Position')
subplot(3,1,3), hold on, grid on, box on
plot(T,v,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(4:5,:))
set(gca,'Xlim',Tlim)
ylabel('Velocity')

figure, clf
subplot(2,1,1), hold on, grid on, box on
plot(T,th,'ko','MarkerFaceColor','k','HandleVisibility','off')
plot(result.t,result.x(6,:))
set(gca,'Xlim',Tlim)
ylabel('Attitude angle')
subplot(2,1,2), hold on, grid on, box on
plot(T,w,'ko','MarkerFaceColor','k','HandleVisibility','off') 
plot(result.t,result.x(7,:))
set(gca,'Xlim',Tlim)
ylabel('Angular rate')

figure, clf
subplot(2,1,1), hold on, grid on, box on
plot(T,U(1,:),'b-o')
plot([0 T(auxdata.kp(1))],[O.bnds.u_min(1) O.bnds.u_min(1)],'r--')
plot([T(auxdata.kp(1)) T(auxdata.kp(2))],3.*[O.auxdata.Tmin O.auxdata.Tmin],'r--')
plot([T(auxdata.kp(2)) T(end)],[O.auxdata.Tmin O.auxdata.Tmin],'r--')

plot([0 T(auxdata.kp(1))],[O.bnds.u_max(1) O.bnds.u_max(1)],'r--')
plot([T(auxdata.kp(1)) T(auxdata.kp(2))],3.*[O.auxdata.Tmax O.auxdata.Tmax],'r--')
plot([T(auxdata.kp(2)) T(end)],[O.auxdata.Tmax O.auxdata.Tmax],'r--')

set(gca,'Xlim',Tlim,'Ylim',[0,ceil(O.bnds.u_max(1))])
ylabel('Thrust')
subplot(2,1,2), hold on, grid on, box on
plot(T,rad2deg(U(2,:)),'b-o')
plot([0 T(end)],rad2deg([O.bnds.u_min(2) O.bnds.u_min(2)]),'r--')
plot([0 T(end)],rad2deg([O.bnds.u_max(2) O.bnds.u_max(2)]),'r--')
set(gca,'Xlim',Tlim,'Ylim',[1.5*(rad2deg(O.bnds.u_min(2))),1.5*(rad2deg(O.bnds.u_max(2)))])
ylabel('Gimbal Angle')
xlabel('Time')

figure, clf, hold on, grid on, box on, axis equal
plot(r(2,:),r(1,:),'ko','MarkerFaceColor','k')
plot(result.x(3,:),result.x(2,:),'b')
ylabel('$z_{\mathcal{I}}$')
xlabel('$y_{\mathcal{I}}$')
title('Landing Trajectory')