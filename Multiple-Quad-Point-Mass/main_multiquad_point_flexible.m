function main_multiquad_point_flexible()
% 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 27  -May-2016
% ====================================================================

%% INITIALZING WORKSPACE
% =====================
% profile on
% Clear workspace
clear; close all; %clc;

% Add Paths
% Geometric Control Toolbox
% addpath('../GeoControl-Toolbox/');

%% INITIALZING PARAMETERS
% ======================
% System constants and parameters
data.nQ = 3; % No. of Quadrotors
data.mL = 0.3; % Mass of the suspended Load
data.g = 9.81 ;
data.e1 = [1;0;0] ;
data.e2 = [0;1;0] ;
data.e3 = [0;0;1] ;

% ---- Quad-1 ----------%
data.params(1).mQ = 0.5 ;
data.params(1).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(1).g = 9.81 ;
data.params(1).e1 = [1;0;0] ;
data.params(1).e2 = [0;1;0] ;
data.params(1).e3 = [0;0;1] ;

data.params(1).m = [0.01,0.01,0.01,0.01,0.5]; % 0.1*ones(1,5); % mass of each link
data.params(1).l =  1; 0.25*ones(1,5); % length of each link 
data.params(1).n = 1; % No. of links suspended

% ---- Quad-2 ----------%
data.params(2).mQ = 0.6 ;
data.params(2).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(2).g = 9.81 ;
data.params(2).e1 = [1;0;0] ;
data.params(2).e2 = [0;1;0] ;
data.params(2).e3 = [0;0;1] ;

data.params(2).m = [0.01,0.01,0.01,0.01,0.5]; % 0.1*ones(1,5); % mass of each link
data.params(2).l = 1;0.25*ones(1,5); % length of each link 
data.params(2).n = 1; % No. of links suspended

% ---- Quad-3 ----------%
data.params(3).mQ = 0.7 ;
data.params(3).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(3).g = 9.81 ;
data.params(3).e1 = [1;0;0] ;
data.params(3).e2 = [0;1;0] ;
data.params(3).e3 = [0;0;1] ;

data.params(3).m = [0.01,0.01,0.01,0.01,0.3]; % 0.1*ones(1,5); % mass of each link
data.params(3).l = 1;0.25*ones(1,5); % length of each link 
data.params(3).n = 1; % No. of links suspended

time = [0:0.1:10]';
for i = 1:length(time)
    fprintf('iter %d in %d\n',i,length(time));
    
    [trajd]= get_desired_traj_multiQuad_pointmass(get_flat_traj(time(i)),data);
    
    for j = 1:data.nQ
        n = data.params(j).n;
        xd_ = [];
        xd_ = [trajd(j).xL(n).x;trajd(j).xL(n).dx{1}];
        xd_ = [xd_;reshape(trajd(j).R, 9,1);trajd(j).Omega];
        for ii = 1:n
            xd_ = [xd_;trajd(j).q(ii).q];
        end
        for ii = 1:n
            xd_ = [xd_;vec_cross(trajd(j).q(ii).q,trajd(j).q(ii).dq{1})];
        end

        xd(i,:,j) = xd_'; 
    
    end
    
end

% ANIMATION
% =========
keyboard;
animate_3dmultiquad_flexcable(time, xd, data.params);


%% INTIALIZING - INTIAL CONDITIONS
% % ================================
% % Zero Position 
% % [trajd0] = get_ini_cond(data.params);
% 
% % Zero Initial Error- Configuration
% % [trajd0] = get_desired_traj(get_flat_traj(0),data.params);
% [trajd0] = get_desired_traj(get_flat_traj_10link(0),data.params);
% 
% xL0 = trajd0.xL(data.params.n).x ;
% vL0 = trajd0.xL(data.params.n).dx{1};
% % qi = {[ 0; sin(th1); cos(th1)],[ 0; sin(th2); cos(th2)] };
% % qi = {[0;0;-1], [0;0;-1], [0;0;-1],[0;0;-1],[0;0;-1] };
% 
% % Make changes to initial conditions Roll-Pitch-Yaw to Rotation Matrix
% % Z-X-Y euler angles
% % R = RPYtoRot_ZXY(0*pi/180, 0*pi/180, 0*pi/180) ;
% R0 = trajd0.R;
% Omega0 = trajd0.Omega;
% 
% q0 = [];
% dq0 = [];
% omega0 = [];
% for i = 1:data.params.n
%     q0 = [q0,trajd0.q(i).q];
%     dq0 = [dq0, trajd0.q(i).dq{1}];
%     omega0 = [omega0, vec_cross(trajd0.q(i).q,trajd0.q(i).dq{1})];
% end
% 
% % state - structure
% % [xL; vL; R; Omega; qi; omegai] - for = 1:n
% % setting up x0 (initial state)
% x0 = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(omega0,numel(omega0),1)];

%% ODE & DIFF-FLAT COMPARISION
% (IGNORE THIS SECTION)
% Times = [0,0.5, 1.78, 5, 8, 10, 15, 100];
% for Time = Times
%     traj_diff = get_desired_traj(get_flat_traj(Time),data.params);
%     xL0 = traj_diff.xL(data.params.n).x ;
%     vL0 = traj_diff.xL(data.params.n).dx{1};
%     R0 = traj_diff.R;
%     Omega0 = traj_diff.Omega;
%     q0 = [];
%     dq0 = [];
%     omega0 = [];
%     for i = 1:data.params.n
%         q0 = [q0,traj_diff.q(i).q];
%         dq0 = [dq0, traj_diff.q(i).dq{1}];
%         omega0 = [omega0, vec_cross(traj_diff.q(i).q,traj_diff.q(i).dq{1})];
%     end
% 
%     dx_diff = [vL0; traj_diff.xL(data.params.n).dx{2}];
%     dx_diff = [dx_diff;  reshape(traj_diff.dR,9,1); traj_diff.dOmega];
%     q_diff = [];
%     dq_diff = [];
%     domega_diff = [];
%     for i = 1:data.params.n
%         q_diff = [q_diff,traj_diff.q(i).q];
%         dq_diff = [dq_diff, traj_diff.q(i).dq{1}];
%         domega_diff = [domega_diff, vec_cross(traj_diff.q(i).q,traj_diff.q(i).dq{2})];
%     end
%     dx_diff = [dx_diff; reshape(dq_diff,numel(dq_diff),1); reshape(domega_diff,numel(domega_diff),1) ];
% 
%     x0 = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(omega0,numel(omega0),1)];
%     [dx_ode, ~, ~,~] = odefun_singleQuadFlex(Time,x0,data);
% 
%     compare_dx = [dx_diff,dx_ode];
%     fprintf('Time = %f, error = %f\n',Time,norm(dx_diff-dx_ode));
%     norm(dx_diff-dx_ode)
%     
%     
% end
% keyboard;

%% SIMULATION
% ==========
% disp('Simulating...') ;
% odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% % odeopts = [] ;
% [t, x] = ode15s(@odefun_singleQuadFlex, [0 1], x0, odeopts, data) ;
% 
% % Computing Various Quantities
% disp('Computing...') ;
% ind = round(linspace(1, length(t), round(1*length(t)))) ;
% % ind = 0:length(t);
% for i = ind
%    [~,xd_] =  odefun_singleQuadFlex(t(i),x(i,:)',data);
%    xd(i,:) = xd_';
%    psi_exL(i) = norm(x(i,1:3)-xd(i,1:3));
%    psi_evL(i) = norm(x(i,4:6)-xd(i,4:6));
% end

%% PLOTS
% % =====
%     figure;
%     subplot(2,2,1);
%     plot(t(ind),x(ind,1),'-g',t(ind),xd(ind,1),':r');
%     grid on; title('x');legend('x','x_d');%axis equal;
%     xlabel('time');ylabel('x [m]');
%     subplot(2,2,2);
%     plot(t(ind),x(ind,2),'-g',t(ind),xd(ind,2),':r');
%     grid on; title('y');legend('y','y_d');%axis equal;
%     xlabel('time');ylabel('y [m]');
%     subplot(2,2,3);
%     plot(t(ind),x(ind,3),'-g',t(ind),xd(ind,3),':r');
%     grid on; title('z');legend('z','z_d');%axis equal;
%     xlabel('time');ylabel('z [m]');
%     subplot(2,2,4);
%     plot3(x(ind,1),x(ind,2),x(ind,3),'-g',xd(ind,1),xd(ind,2),xd(ind,3),':r');
%     grid on; title('trajectory');legend('traj','traj_d');%axis equal;
%     xlabel('x-axis');ylabel('y-axis');zlabel('z-axis');
% 
%     figure;
%     subplot(2,1,1);
%     plot(t(ind),psi_exL(ind));
%     grid on; title('position error');legend('psi-exL');
%     subplot(2,1,2);
%     plot(t(ind),psi_evL(ind));
%     grid on; title('velocity error');legend('psi-evL');
% 
% % ANIMATION
% % =========
% keyboard;
% animate_3dquad_flexcable(t(ind), xd(ind,:), data.params);
% 
% profsave;

end

%%
function[dx, xd, f,M] = odefun_singleQuadFlex(t,x,data)
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
% [trajd] = get_desired_traj(get_flat_traj(t),data.params);
[trajd] = get_desired_traj(get_flat_traj_10link(t),data.params);
% inputs
f = trajd.f;
M = trajd.M;
% f = M00*g;
% M = zeros(3,1);

% Extracing states
% ----------------
xL = x(1:3);
vL = x(4:6);
R = reshape(x(7:15),3,3);
Omega = x(16:18);
q = reshape(x(19:3*n+18),3,n);
omega = reshape(x(3*n+19:6*n+18),3,n);


    % Equations of Motion
    % -------------------
    dx = [];

    xL_dot = vL;
    q_dot = zeros(size(q));
    for i = 1:n
       q_dot(:,i) = vec_cross(omega(:,i),q(:,i));
    end

%     lhsMat = [];
    lhsMat = M00*eye(3);
%     for ii = 1:n
%         lhsMat = [lhsMat, -M0i(ii)*hat_map(q(:,ii))];
%     end
    tmp_l = zeros(3*n,3*n);
    tmp_l2 = [];
    for ii = 1:n
        lhsMat = [lhsMat, -M0i(ii)*hat_map(q(:,ii))];
        tmp_l2 = [tmp_l2; hat_map(q(:,ii))*Mi0(ii)];
        for jj = 1:n
            if ii==jj
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = Mij(ii,jj)*eye(3);
            else
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = -Mij(ii,jj)*hat_map(q(:,ii))*hat_map(q(:,jj));
            end
        end
    end
    
%     for i = 1:n
%      tmp_l2 = [tmp_l2; hat_map(q(:,i))*Mi0(i)];   
%     end
    lhsMat = [lhsMat;tmp_l2, tmp_l];
    
%     rhsMat = [];
    rhsMat = f*R*e3 - M00*g*e3;
    for i = 1:n
       rhsMat = rhsMat + M0i(i)*norm(omega(:,i))^2*q(:,i);  
    end
    for i = 1:n
        tmp_rhs =  - sum(m(i:n))*g*l(i)*hat_map(q(:,i))*e3;
        for j = 1:n
            if j~=i
               tmp_rhs = tmp_rhs + Mij(i,j)*norm(omega(:,j))^2*hat_map(q(:,i))*q(:,j) ;
            end
        end
        rhsMat = [rhsMat; tmp_rhs];
    end
    
    d2x = lhsMat\rhsMat;
    omega_dot = reshape(d2x(4:end),3,n);
    
        
    vL_dot = d2x(1:3,1);
    for i = 1:n
        vL_dot = vL_dot + l(i)*(vec_cross(omega_dot(:,i),q(:,i))+ vec_cross(omega(:,i),q_dot(:,i)));
    end
    
    % Quadrotor Attitude
    R_dot = R*hat_map(Omega) ;
    Omega_dot = (J)\( -vec_cross(Omega, J*Omega) + M ) ;
    
% Computing xd
% ------------
xd = [trajd.xL(n).x;trajd.xL(n).dx{1}];
xd = [xd;reshape(trajd.R, 9,1);trajd.Omega];
for i = 1:n
    xd = [xd;trajd.q(i).q];
end
for i = 1:n
    xd = [xd;vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
end
% Computing dx
%-------------
dx = [xL_dot;
      vL_dot;
      reshape(R_dot, 9,1) ;
      Omega_dot;
      reshape(q_dot,numel(q_dot),1);
      reshape(omega_dot,numel(omega_dot),1)];

if nargout <= 1
   fprintf('Sim time %0.4f seconds \n',t);
end
    
end


%%
function[trajd] = get_ini_cond(params)

traj.x = zeros(3,1);

for i = 1:(4+2*params.n)
    traj.dx{i} = zeros(3,1);
end
   
[trajd]= get_desired_traj(traj,params);
th = 178*pi/180;
trajd.q(params.n).q = [0;sin(th);cos(th)];


end
