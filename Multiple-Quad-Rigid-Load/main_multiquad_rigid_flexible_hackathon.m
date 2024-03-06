function main_multiquad_rigid_flexible_hackathon()
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
data.vert = getVertices(2,2,0.5)-[0,0,0.25];
data.r = data.vert(5:end,:)';%[RotZ(0)*[1;0;0], RotZ(-30*pi/180)*[2;0;0], RotZ(-165*pi/180)*[1.5;0;0], RotZ(45*pi/180)*[2.5;0;0]];
data.RL = RotX(5*pi/180);%[1,0,0;0,1,0;0,0,1];
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
data.params(1).l = 1.25*ones(1,5); % length of each link 
data.params(1).n = 1; % No. of links suspended

% ---- Quad-2 ----------%
data.params(2).mQ = 0.5 ;
data.params(2).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(2).g = 9.81 ;
data.params(2).e1 = [1;0;0] ;
data.params(2).e2 = [0;1;0] ;
data.params(2).e3 = [0;0;1] ;

data.params(2).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(2).l = 1.25*ones(1,5); % length of each link 
data.params(2).n = 1; % No. of links suspended

% ---- Quad-3 ----------%
data.params(3).mQ = 0.5 ;
data.params(3).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(3).g = 9.81 ;
data.params(3).e1 = [1;0;0] ;
data.params(3).e2 = [0;1;0] ;
data.params(3).e3 = [0;0;1] ;

data.params(3).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(3).l = 1.25*ones(1,5); % length of each link 
data.params(3).n = 1; % No. of links suspended

% ---- Quad-3 ----------%
data.params(4).mQ = 0.5 ;
data.params(4).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(4).g = 1.81 ;
data.params(4).e1 = [1;0;0] ;
data.params(4).e2 = [0;1;0] ;
data.params(4).e3 = [0;0;1] ;

data.params(4).m = [0.01,0.01,0.01,0.01,0.0]; % 0.1*ones(1,5); % mass of each link
data.params(4).l = 1.25*ones(1,5); % length of each link 
data.params(4).n = 1; % No. of links suspended

% ----------------------%
data.MT  = data.mL;
for i = 1:data.nQ
    data.MT = data.MT + sum(data.params(i).m(1:data.params(i).n));
end
data.MiT = @(i) data.params(i).mQ + sum(data.params(i).m(1:data.params(i).n));
data.M0ij = @(i,j) data.params(i).mQ + (j>1)*sum(data.params(i).m(1:(j-1)));
data.Mij = @(i,j,k) sum(data.params(i).m(max(j,k):data.params(i).n));

data.lij = [];
data.mij = [];
data.ni = [];

for i = 1:data.nQ
    data.lij = [data.lij; data.params(i).l];
    data.mij = [data.mij; data.params(i).m];
    data.ni =  [data.ni; data.params(i).n]; 
end

%% ANIMATION - from DIFF-FLAT TRAJECTORY 
% =========================================================================
[trajd]= get_desired_traj_multiQuad_rigidload(get_flat_traj(0),data);

time = [0:0.1:1]';
for i = 1:length(time)
    fprintf('iter %d in %d\n',i,length(time));
    
    [trajd]= get_desired_traj_multiQuad_rigidload(get_flat_traj(time(i)),data);
    
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
    
        xLcom_ = trajd(j).xLcom.x;
        xLcom(i,:,j) = xLcom_';
        
    end
    
end

% ANIMATION
% =========
keyboard;
animate_3dmultiquad_flexcable_rigid_hackathon(time, xd, xLcom, data);


%% INTIALIZING - INTIAL CONDITIONS
% % ================================
% Zero Position 
% [trajd0] = get_ini_cond(data.params);

% Zero Initial Error- Configuration
% [trajd0] = get_desired_traj(get_flat_traj(0),data.params);
% [trajd0] = get_desired_traj_multiQuad_rigidload(get_flat_traj(0),data);
% 
% xL0 = trajd0(1).xLcom.x ;
% vL0 = trajd0(1).xLcom.dx{1};
% RL0 = trajd0(1).RL;
% OmegaL0 = trajd0(1).OmegaL;
% % qi = {[ 0; sin(th1); cos(th1)],[ 0; sin(th2); cos(th2)] };
% % qi = {[0;0;-1], [0;0;-1], [0;0;-1],[0;0;-1],[0;0;-1] };
% 
% % Make changes to initial conditions Roll-Pitch-Yaw to Rotation Matrix
% % Z-X-Y euler angles
% % R = RPYtoRot_ZXY(0*pi/180, 0*pi/180, 0*pi/180) ;
% 
% % state - structure
% % [xL; vL; RL; OmegaL; Ri, Omegai qij; omegaij] - for i = 1:nQ & j = 1:ni
% % setting up x0 (initial state)
% x0 = [xL0; vL0; reshape(RL0,9,1); OmegaL0]; 
% 
% R0 = [];
% Omega0 = [];
% q0 = [];
% dq0 = [];
% omega0 = [];
% state_q = [];
% 
% for i = 1:data.nQ
%     R0 = trajd0(i).R;
%     Omega0= trajd0(i).Omega;
%     x0 = [x0; reshape(R0,9,1); Omega0];
%     
%     q0 = [];
%     dq0 = [];
%     omega0 = [];
%     for j = 1:data.params(i).n
%         q0 = [q0,trajd0(i).q(j).q];
%         dq0 = [dq0, trajd0(i).q(j).dq{1}];
%         omega0 = [omega0, vec_cross(trajd0(i).q(j).q,trajd0(i).q(j).dq{1})];
%     end
%     state_q = [state_q; reshape(q0, numel(q0),1); reshape(omega0,numel(omega0),1)];
% end
% 
% x0 = [x0;state_q];



%% SIMULATION
% ==========
% disp('Simulating...') ;
% odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% % odeopts = [] ;
% [t, x] = ode15s(@odefun_multiQuad_rigid, [0 .1], x0, odeopts, data) ;
% % 
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
function[dx, xd, f,M] = odefun_multiQuad_rigid(t,x,data)
% Extracing parameters
% --------------------
% System constants and parameters
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1;
e2 = data.e2;
e3 = data.e3 ;
r = data.r;
RL = data.RL;
JL = data.JL;

MT = data.MT;
MiT = data.MiT;
M0ij = data.M0ij;
Mij = data.Mij;

mij = data.mij;
lij = data.lij;
ni = data.ni;

% fetching desired states
% -----------------------
% [trajd] = get_desired_traj(get_flat_traj(t),data.params);
[trajd] = get_desired_traj_multiQuad_rigidload(get_flat_traj(t),data);


% state - structure
% [xL; vL; RL; OmegaL; (Ri; Omegai); (qij; omegaij)] - for i = 1:nQ & j = 1:ni

% Extracing states
% ----------------
xL = x(1:3);
vL = x(4:6);
RL = reshape(x(7:15),3,3);
OmegaL = x(16:18);

tempR = x(19:(12*nQ+18));
tempq = x((12*nQ+19):end);

for i = 1:nQ
    Ri(:,:,i) = reshape(tempR(12*i-11:12*i-3),3,3);
    Omegai(:,1,i) = tempR(12*i-2:12*i);
    
    tempq2 = tempq((1-6*ni(i)+6*sum(ni(1:i))): 6*sum(ni(1:i)));
    
    qij{i} = reshape(tempq2(1:3*ni(i)),3,ni(i));
    omegaij{i} = reshape(tempq2(1+3*ni(i):6*ni(i)),3,ni(i));
    
    fi(1,i) = trajd(i).f;
    Mi(:,i) = trajd(i).M;
    
end


    % Equations of Motion
    % ===================
    dx = [];
    
    N_x0_Omega0 = zeros(3);
    
    for i = 1:nQ
        N_x0_Omega0 = N_x0_Omega0 - MiT(i)*RL*hat_map(r(:,i));

        N_x0i{i} = [];
        N_Omega0i{i} = [];
        N_ix0{i} = [];
        N_iOmega0{i} = [];
        N_qqi{i} = zeros(ni(i)*3);
        
        for j = 1:ni(i)
           N_x0i{i} = [N_x0i{i}, -M0ij(i,j)*lij(i,j)*eye(3)];
           N_Omega0i{i} = [N_Omega0i{i}, -M0ij(i,j)*lij(i,j)*hat_map(r(:,i))*RL'];
           N_ix0{i} = [N_ix0{i}, -M0ij(i,j)*hat_map(qij{i}(:,j))^2];
           N_iOmega0{i} = [N_iOmega0{i}, M0ij(i,j)*hat_map(qij{i}(:,j))^2*RL*hat_map(r(:,i))];
           
           for k = 1:ni(i)
               if j == k
                   N_qqi{i}(3*j-2:3*j,3*k-2:3*k) = -Mij(i,j,k)*lij(i,k)*eye(3);
               else
                   N_qqi{i}(3*j-2:3*j,3*k-2:3*k) = Mij(i,j,k)*lij(i,k)*hat_map(qij{i}(:,k))^2;
               end
           end
        end    
        
        N_ix0{i} = transpose(N_ix0{i});
        N_iOmega0{i} = transpose(N_iOmega0{i});
        
    end
    
    N_Omega0_x0 = N_x0_Omega0';
 
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
