function verify_diff_flat()
%
% Date: June-10-2016
% Last Updated: June-10-2016
%
%% INITIALZING WORKSPACE
% =====================
% profile on
% Clear workspace
clear; close all; %clc;

% Add Paths
% Geometric Control Toolbox
addpath('../GeoControl-Toolbox/');

%% INITIALZING PARAMETERS
% ======================
% System constants and parameters
data.nQ = 4; % No. of Quadrotors
data.mL = 0.3; % Mass of the suspended Load
data.g = 9.81 ;
data.e1 = [1;0;0] ;
data.e2 = [0;1;0] ;
data.e3 = [0;0;1] ;
% data.r = [[0.5;0;0],[0.1732;-0.1;0],[-0.6761;-0.1812;0],[0.2121;0.2121;0]];
data.r = [RotZ(0)*[1;0;0], RotZ(-30*pi/180)*[2;0;0], RotZ(-165*pi/180)*[1.5;0;0], RotZ(45*pi/180)*[2.5;0;0]];
data.RL = eye(3);%RotX(5*pi/180);%[1,0,0;0,1,0;0,0,1];
data.JL = diag([2.1, 1.87, 3.97]*10e-2);
data.Phi = [];
for i = 1:data.nQ
    data.Phi = [data.Phi,[eye(3);hat_map(data.r(:,i))]];
end
data.Phi_pinv = pinv(data.Phi);
data.N = [];
r = [data.r, data.r(:,1)];
u = cell(1,1);
for i = 1:(data.nQ)
   for j = 1:data.nQ
       if i~=j
           u{i,j} = (r(:,j)-r(:,i))/norm(r(:,j)-r(:,i));
       end
   end
end
% (FIX THIS - write code to automate following variable - data.N)
% ----------------------------------------------------------------
data.N = [[u{1,2};-u{1,2};zeros(3,1);zeros(3,1)],...
            [u{1,3};zeros(3,1);-u{1,3};zeros(3,1)],...
            [u{1,4};zeros(3,1);zeros(3,1);-u{1,4}],...
            [zeros(3,1);u{2,3};-u{2,3};zeros(3,1)],...
            [zeros(3,1);u{2,4};zeros(3,1);-u{2,4}],...
            [zeros(3,1);zeros(3,1);u{3,4};-u{3,4}]];

% ---- Quad-1 ----------%
data.params(1).mQ = 0.5 ;
data.params(1).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(1).g = 9.81 ;
data.params(1).e1 = [1;0;0] ;
data.params(1).e2 = [0;1;0] ;
data.params(1).e3 = [0;0;1] ;

data.params(1).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(1).l = 0.25*ones(1,5); % length of each link 
data.params(1).n = 5; % No. of links suspended

% data.params(1).M00 = data.params(1).mQ + sum(data.params(1).m(1:data.params(1).n));
% data.params(1).M0i = @(i) sum(data.params(1).m(i:data.params(1).n))*data.params(1).l(i); 
% data.params(1).Mi0 = data.params(1).M0i;
% data.params(1).Mij = @(i,j) sum(data.params(1).m(max(i,j):data.params(1).n))*data.params(1).l(i)*data.params(1).l(j);

% ---- Quad-2 ----------%
data.params(2).mQ = 0.5 ;
data.params(2).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(2).g = 9.81 ;
data.params(2).e1 = [1;0;0] ;
data.params(2).e2 = [0;1;0] ;
data.params(2).e3 = [0;0;1] ;

data.params(2).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(2).l = 0.25*ones(1,5); % length of each link 
data.params(2).n = 5; % No. of links suspended

% data.params(2).M00 = data.params(2).mQ + sum(data.params(2).m(1:data.params(2).n));
% data.params(2).M0i = @(i) sum(data.params(2).m(i:data.params(2).n))*data.params(2).l(i); 
% data.params(2).Mi0 = data.params(2).M0i;
% data.params(2).Mij = @(i,j) sum(data.params(2).m(max(i,j):data.params(2).n))*data.params(2).l(i)*data.params(2).l(j);

% ---- Quad-3 ----------%
data.params(3).mQ = 0.5 ;
data.params(3).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(3).g = 9.81 ;
data.params(3).e1 = [1;0;0] ;
data.params(3).e2 = [0;1;0] ;
data.params(3).e3 = [0;0;1] ;

data.params(3).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(3).l = 0.25*ones(1,5); % length of each link 
data.params(3).n = 5; % No. of links suspended

% data.params(3).M00 = data.params(3).mQ + sum(data.params(3).m(1:data.params(3).n));
% data.params(3).M0i = @(i) sum(data.params(3).m(i:data.params(3).n))*data.params(3).l(i); 
% data.params(3).Mi0 = data.params(3).M0i;
% data.params(3).Mij = @(i,j) sum(data.params(3).m(max(i,j):data.params(3).n))*data.params(3).l(i)*data.params(3).l(j);

% ---- Quad-3 ----------%
data.params(4).mQ = 0.5 ;
data.params(4).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(4).g = 9.81 ;
data.params(4).e1 = [1;0;0] ;
data.params(4).e2 = [0;1;0] ;
data.params(4).e3 = [0;0;1] ;

data.params(4).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(4).l = 0.25*ones(1,5); % length of each link 
data.params(4).n = 5; % No. of links suspended

% data.params(4).M00 = data.params(4).mQ + sum(data.params(4).m(1:data.params(4).n));
% data.params(4).M0i = @(i) sum(data.params(4).m(i:data.params(4).n))*data.params(4).l(i); 
% data.params(4).Mi0 = data.params(4).M0i;
% data.params(4).Mij = @(i,j) sum(data.params(4).m(max(i,j):data.params(4).n))*data.params(4).l(i)*data.params(4).l(j);


%% load
% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% xQ0 = zeros(3,1);
% vQ0 = zeros(3,1);
% 
% R0 = RPYtoRot_ZXY(0*pi/180,0*pi/180, 0*pi/180) ;
% Omega0 = zeros(3,1);

% Zero Initial Error- Configuration
% ---------------------------------
[trajd0] = get_desired_traj_multiQuad_rigidload(get_flat_traj(0),data);
xQ0 = trajd0(1).xLcom.x;
vQ0 = trajd0(1).xLcom.dx{1};

R0 = trajd0(1).RL;
Omega0 = trajd0(1).OmegaL;

% state - structure
% -----------------
% [xL; vL; R; Omega]
% setting up x0 (initial state)
% -----------------------------
x0 = [xQ0; vQ0; reshape(R0,9,1); Omega0 ];


% SIMULATION
% ==========
disp('Simulating...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% odeopts = [] ;
[t, x] = ode15s(@odefun_load_dynamics, [0 25], x0, odeopts, data) ;
% 
% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(1*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_] =  odefun_load_dynamics(t(i),x(i,:)',data);
   xd(i,:) = xd_';
   psi_exL(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evL(i) = norm(x(i,4:6)-xd(i,4:6));
   psi_R(i) = 0.5*trace(eye(3,3) - reshape(xd_(7:15),3,3)'*reshape(x(i,7:15),3,3));
end


% PLOTS
% =====
    figure;
    subplot(2,2,1);
    plot(t(ind),x(ind,1),'-g',t(ind),xd(ind,1),':r');
    grid on; title('x');legend('x','x_d');%axis equal;
    xlabel('time');ylabel('x [m]');
    subplot(2,2,2);
    plot(t(ind),x(ind,2),'-g',t(ind),xd(ind,2),':r');
    grid on; title('y');legend('y','y_d');%axis equal;
    xlabel('time');ylabel('y [m]');
    subplot(2,2,3);
    plot(t(ind),x(ind,3),'-g',t(ind),xd(ind,3),':r');
    grid on; title('z');legend('z','z_d');%axis equal;
    xlabel('time');ylabel('z [m]');
    subplot(2,2,4);
    plot3(x(ind,1),x(ind,2),x(ind,3),'-g',xd(ind,1),xd(ind,2),xd(ind,3),':r');
    grid on; title('trajectory');legend('traj','traj_d');%axis equal;
    xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');

    figure;
    subplot(2,2,1);
    plot(t(ind),psi_exL(ind));
    grid on; title('position error');legend('psi-exL');
    subplot(2,2,2);
    plot(t(ind),psi_evL(ind));
    grid on; title('velocity error');legend('psi-evL');
    subplot(2,2,3);
    plot(t(ind),psi_R(ind));
    grid on; title('Rotation error');legend('psi-R');
    
%% quadrotor
qNo = 1;
data.qNo = qNo;
% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% [trajd0] = get_ini_cond(data.params);

% Zero Initial Error- Configuration
% ---------------------------------
[trajd0] = get_desired_traj_multiQuad_rigidload(get_flat_traj(0),data);
xL0 = trajd0(qNo).xL(data.params(qNo).n-1).x ;
vL0 = trajd0(qNo).xL(data.params(qNo).n-1).dx{1};

R0 = trajd0(qNo).R;
Omega0 = trajd0(qNo).Omega;

q0 = [];
dq0 = [];
omega0 = [];
for i = 1:(data.params(qNo).n-1)
    q0 = [q0,trajd0(qNo).q(i).q];
    dq0 = [dq0, trajd0(qNo).q(i).dq{1}];
    omega0 = [omega0, vec_cross(trajd0(qNo).q(i).q,trajd0(qNo).q(i).dq{1})];
end

% state - structure
% -----------------
% [xL; vL; R; Omega; qi; omegai] - for = 1:n
% setting up x0 (initial state)
% -----------------------------
x0 = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(dq0,numel(dq0),1)];

% SIMULATION
% ==========
disp('Simulating...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% odeopts = [] ;
[t, x] = ode15s(@odefun_quad_dynamics, [0 15], x0, odeopts, data) ;

% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(1*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_] =  odefun_quad_dynamics(t(i),x(i,:)',data);
   xd(i,:) = xd_';
   psi_exL(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evL(i) = norm(x(i,4:6)-xd(i,4:6));
   psi_R(i) = 0.5*trace(eye(3,3) - reshape(xd_(7:15),3,3)'*reshape(x(i,7:15),3,3));
end

% PLOTS
% =====

    figure;
    subplot(2,2,1);
    plot(t(ind),psi_exL(ind));
    grid on; title('position error');legend('psi-exL');
    subplot(2,2,2);
    plot(t(ind),psi_evL(ind));
    grid on; title('velocity error');legend('psi-evL');
    subplot(2,2,3);
    plot(t(ind),psi_R(ind));
    grid on; title('Rotation error');legend('psi-R');
    subplot(2,2,4);
    plot3(x(ind,1),x(ind,2),x(ind,3),'-g',xd(ind,1),xd(ind,2),xd(ind,3),':r');
    grid on; title('trajectory');legend('traj','traj_d');%axis equal;
    xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');

end

function[dx,xd] = odefun_load_dynamics(t,x,data)
% Parameters
% ----------
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1;
e2 = data.e2;
e3 = data.e3;
r = data.r;
JL = data.JL;

% Fetching desired Trajectory
% ---------------------------
[trajd]= get_desired_traj_multiQuad_rigidload(get_flat_traj(t),data);



% Extracting state
% ----------------
xL = x(1:3);
vL = x(4:6);
RL = reshape(x(7:15),3,3);
OmegaL = x(16:18);


sum_Tensions = zeros(3,1);
sum_Torques = zeros(3,1);
for i = 1:length(trajd)
    sum_Tensions = sum_Tensions + trajd(i).B(end).B;
    sum_Torques = sum_Torques - vec_cross(r(:,i),RL'*trajd(i).B(end).B);
end

    % Equations of Motion
    % -------------------
    xL_dot = vL;
    RL_dot = RL*hat_map(OmegaL);
    
    vL_dot = sum_Tensions/mL - g*e3; 
    OmegaL_dot = JL\(-vec_cross(OmegaL,JL*OmegaL)+sum_Torques);
    
    
    
% Computing xd
% -----------
xd = [trajd(1).xLcom.x;
        trajd(1).xLcom.dx{1};
        reshape(trajd(1).RL,9,1);
        trajd(1).OmegaL];
    
% Computing dx
% ------------
dx = [xL_dot;
        vL_dot;
        reshape(RL_dot,9,1);
        OmegaL_dot];


if nargout <= 1
   fprintf('Sim time %0.4f seconds \n',t);
end


end


function[dx,xd] = odefun_quad_dynamics(t,x,data)
% this function uses 'q & dq' as states instead of 'q & omega' which was
% the case in 'single-Quad-Flex-Cable>main_quad_flexible>odefun_singleQuadFlex'
% Extracing parameters
% --------------------
qNo = data.qNo;
% Dynamics of quadrotor suspended with load Constants
mQ = data.params(qNo).mQ;
J = data.params(qNo).J;
g = data.params(qNo).g;
e1 = data.params(qNo).e1;
e2 = data.params(qNo).e2;
e3 = data.params(qNo).e3;

m = data.params(qNo).m; % mass of each link
l = data.params(qNo).l; % length of each link 
n = data.params(qNo).n-1; % No. of links suspended

M00 = mQ + sum(m(1:n));
M0i = @(i) sum(m(i:n))*l(i); 
Mi0 = M0i;
Mij = @(i,j) sum(m(max(i,j):n))*l(i)*l(j);

% fetching desired states
% -----------------------
[trajd]= get_desired_traj_multiQuad_rigidload(get_flat_traj(t),data);
f = trajd(qNo).f;
M = trajd(qNo).M;
T = -trajd(qNo).B(end).B;

% Extracing states
% ----------------
xL = x(1:3);
vL = x(4:6);
R = reshape(x(7:15),3,3);
Omega = x(16:18);
q = reshape(x(19:3*n+18),3,n);
dq = reshape(x(3*n+19:6*n+18),3,n);


    % Equations of Motion
    % -------------------
    dx = [];

    xL_dot = vL;
    q_dot = zeros(size(q));
    for i = 1:n
       q_dot(:,i) = dq(:,i);%vec_cross(omega(:,i),q(:,i));
    end

%     lhsMat = [];
    lhsMat = M00*eye(3);
%     for ii = 1:n
%         lhsMat = [lhsMat, -M0i(ii)*hat_map(q(:,ii))];
%     end
    tmp_l = zeros(3*n,3*n);
    tmp_l2 = [];
    for ii = 1:n
        lhsMat = [lhsMat, M0i(ii)*eye(3)];
        tmp_l2 = [tmp_l2; -hat_map(q(:,ii))^2*Mi0(ii)];
        for jj = 1:n
            if ii==jj
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = Mij(ii,jj)*eye(3);
            else
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = -Mij(ii,jj)*hat_map(q(:,ii))^2;
            end
        end
    end
    
%     for i = 1:n
%      tmp_l2 = [tmp_l2; hat_map(q(:,i))*Mi0(i)];   
%     end
    lhsMat = [lhsMat;tmp_l2, tmp_l];
    
%     rhsMat = [];
    rhsMat = f*R*e3 - M00*g*e3 + T;
    for i = 1:n
        tmp_rhs = -norm(dq(:,i))^2*Mij(i,i)*q(:,i) + sum(m(i:n))*g*l(i)*hat_map(q(:,i))^2*e3 - l(i)*hat_map(q(:,i))^2*T;
        rhsMat = [rhsMat; tmp_rhs];
    end
    
    d2x = lhsMat\rhsMat;
    dq_dot = reshape(d2x(4:end),3,n);
    
        
    vL_dot = d2x(1:3,1);
    for i = 1:n
        vL_dot = vL_dot + l(i)*dq_dot(:,i);
    end
    
    % Quadrotor Attitude
    R_dot = R*hat_map(Omega) ;
    Omega_dot = (J)\( -vec_cross(Omega, J*Omega) + M ) ;
    
% Computing xd
% ------------
xd = [ trajd(qNo).xL(n).x;trajd(qNo).xL(n).dx{1}];
xd = [xd;reshape(trajd(qNo).R, 9,1);trajd(qNo).Omega];
for i = 1:n
    xd = [xd;trajd(qNo).q(i).q];
end
for i = 1:n
    xd = [xd;trajd(qNo).q(i).dq{1}];
end
% Computing dx
%-------------
dx = [xL_dot;
      vL_dot;
      reshape(R_dot, 9,1) ;
      Omega_dot;
      reshape(q_dot,numel(q_dot),1);
      reshape(dq_dot,numel(dq_dot),1)];

if nargout <= 1
   fprintf('Sim time %0.4f seconds \n',t);
end
    
end
