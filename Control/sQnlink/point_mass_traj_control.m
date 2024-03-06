function point_mass_traj_control()
%% point_mass_traj_control
% by vkotaru@andrew.cm.edu
% Date: July-5-2016
% Last Updated: July-5-2016
% -----------------------------------------------------------------------
% control of 'quadrotor with suspended flexible cable of 2 links'
% with desired trajectory planned using 'Differential Flatness'
% -----------------------------------------------------------------------
% quad_plan_traj() 
%
% Input Argument: 
% 				None 
% Return:
% 		None
% --

%% INITIALZING WORKSPACE
% =====================
% profile on
% Clear workspace
% clear; 
close all; 
% clc;

% Add Paths
% Geometric Control Toolbox
addpath('../GeoControl-Toolbox/');


%% INITIALZING PARAMETERS
% ======================
% System constants and parameters
data.params.mQ = 0.5 ;
data.params.J = diag([0.557, 0.557, 1.05]*10e-2);
data.params.g = 9.81 ;
data.params.e1 = [1;0;0] ;
data.params.e2 = [0;1;0] ;
data.params.e3 = [0;0;1] ;

data.params.m = [0.5,0.5]; % mass of each link
data.params.l = 0.625*ones(1,2); % length of each link 
data.params.n = 2; % No. of links suspended

% data.params.m = [0.14]; % 0.1*ones(1,5); % mass of each link
% data.params.l = 1.25*ones(1,1); % length of each link 
% data.params.n = 1; % No. of links suspended

n = data.params.n;
data.params.M00 = data.params.mQ + sum(data.params.m(1:data.params.n));
data.params.M0i = @(i) sum(data.params.m(i:data.params.n))*data.params.l(i); 
data.params.Mi0 = data.params.M0i;
data.params.Mij = @(i,j) sum(data.params.m(max(i,j):data.params.n))*data.params.l(i)*data.params.l(j);

%% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% [trajd0] = get_ini_cond(data.params);

% Zero Initial Error- Configuration
% ---------------------------------
[trajd0] = Flat2state(get_flat_traj_10link(0,1),data.params);

xL0 = trajd0.xL(n).x ;zeros(3,1);% trajd0.x ;
vL0 = trajd0.xL(n).dx{1};zeros(3,1);%trajd0.dx{1};
q0 = trajd0.q(n).q;
dq0 = trajd0.q(n).dq{1};

x0 = [xL0;vL0;q0;dq0];


%% SIMULATION
% ==========
disp('Simulating...') ;
odeopts = odeset('RelTol', 1e-12, 'AbsTol', 1e-13) ;
% odeopts = [] ;
[t, x] = ode15s(@odefun_singleQuadFlex, [0 1], x0, odeopts, data) ;

% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(.5*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_] =  odefun_singleQuadFlex(t(i),x(i,:)',data);
   xd(i,:) = xd_';
   psi_exL(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evL(i) = norm(x(i,4:6)-xd(i,4:6));
   
end

%% PLOTS
% ======
figure;
subplot(1,2,1);
plot(t(ind),psi_exL(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Position Error');
legend('psi-exL');
subplot(1,2,2);
plot(t(ind),psi_evL(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Velocity Error');
legend('psi-evL');


figure;
subplot(2,2,1);
plot(t(ind),xd(ind,1),':r',t(ind),x(ind,1),'b'); grid on;
subplot(2,2,2);
plot(t(ind),xd(ind,2),':r',t(ind),x(ind,2),'b'); grid on;
subplot(2,2,3);
plot(t(ind),xd(ind,3),':r',t(ind),x(ind,3),'b'); grid on;
subplot(2,2,4);
plot3(xd(ind,1),xd(ind,2),xd(ind,3),':r',x(ind,1),x(ind,2),x(ind,3),'b'); grid on;

keyboard;
end



%%
function[dx, xd] = odefun_singleQuadFlex(t,x,data)
% Extracing parameters
% --------------------
% Dynamics of quadrotor suspended with load Constants
mQ = data.params.mQ;
J = data.params.J;
g = data.params.g;
e1 = data.params.e1;
e2 = data.params.e2;
e3 = data.params.e3;

m = data.params.m; % mass of each link
l = data.params.l; % length of each link 
n = data.params.n; % No. of links suspended

M00 = data.params.M00;
M0i = data.params.M0i; 
Mi0 = data.params.Mi0;
Mij = data.params.Mij;

% fetching desired states
% -----------------------
% [trajd] = Flat2state(get_flat_traj(t),data.params);
[trajd] =  Flat2state(get_flat_traj_10link(t,1),data.params);


% Extracing states
% ----------------
xL = x(1:3);
vL = x(4:6);
q = x(7:9);
dq = x(10:12);

  	% Control
	  % -------
    eL2 = xL - trajd.xL(n).x;
    deL2 = vL - trajd.xL(n).dx{1};
    epsilon_bar = 0.08 ;
    kp_xy = 0.3/epsilon_bar^2 ; kd_xy = 0.6/epsilon_bar ;
    k1 = diag([kp_xy kp_xy 2]) ; k2 = diag([kd_xy kd_xy 1.5]) ;
    Fpd_2 = -k1*eL2 - k2*deL2;
    Td2 = (Fpd_2 + m(n)*(g*e3 + trajd.xL(n).dx{2}));
    qd2 = -Td2/norm(Td2);
    
    eL1 = (xL -l(n)*q)- (trajd.xL(n).x-l(n)*qd2);
    deL1 = (vL -l(n)*q)- trajd.xL(1).dx{1};
    epsilon_bar = 0.08 ;
    kp_xy = 0.3/epsilon_bar^2 ; kd_xy = 0.6/epsilon_bar ;
    k1 = diag([kp_xy kp_xy 2]) ; k2 = diag([kd_xy kd_xy 1.5]) ;
    Fpd_1 = -k1*eL1 - k2*deL1;
    F = (Fpd_1 + m(1)*(g*e3 + trajd.xL(1).dx{2})) -Td2;
        

    % Equations of Motion
    % -------------------
    dx = [];

    xL_dot = vL;
    vL_dot = (vec_dot(q,F)-m(1)*l(n)*vec_dot(dq,dq))*q/(m(1)+m(2)) -g*e3;
    q_dot = dq;
    
    dq_dot = (1/(m(1)*l(n)))*(vec_cross(q,vec_cross(q,F))) -vec_dot(dq,dq)*q;
    
    
% Computing xd
% ------------
xd = [trajd.xL(n).x;trajd.xL(n).dx{1}];
xd = [xd;trajd.q(n).q; trajd.q(n).dq{1}];
    
% Computing dx
%-------------
dx = [xL_dot;
        vL_dot;
        q_dot;
        dq_dot];

if nargout <= 1
   fprintf('Sim time %0.4f seconds \n',t);
end
    
end

%%

