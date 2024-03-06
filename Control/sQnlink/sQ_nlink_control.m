function sQ_nlink_control()
%% sQ_nlink_control
% by vkotaru@andrew.cm.edu
% Date: July-09-2016
% Last Updated: July-09-2016
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
% addpath(genpath('/home/venkata/git/Lab/Quad_FlexCable/GeoControl-Toolbox'));
addpath(genpath('../../GeoControl-Toolbox'));


%% INITIALZING PARAMETERS
% ======================
% System constants and parameters
data.params.mQ = 0.85 ;
data.params.J = diag([0.557, 0.557, 1.05]*10e-2);
data.params.g = 9.81 ;
data.params.e1 = [1;0;0] ;
data.params.e2 = [0;1;0] ;
data.params.e3 = [0;0;1] ;

data.params.m = [0.1,0.5,0.1,0.1,0.5]; % mass of each link
data.params.l = 0.25*ones(1,5); % length of each link 
data.params.n = 5; % No. of links suspended

% data.params.m = [0.5, 0.5]; % mass of each link
% data.params.l = 0.5*ones(1,2); % length of each link 
% data.params.n = 2; % 

% data.params.m = [0.5]; % mass of each link
% data.params.l = 0.5; % length of each link 
% data.params.n = 1; % 

n = data.params.n;
data.params.M00 = data.params.mQ + sum(data.params.m(1:data.params.n));
data.params.M0i = @(i) sum(data.params.m(i:data.params.n))*data.params.l(i); 
data.params.Mi0 = data.params.M0i;
data.params.Mij = @(i,j) sum(data.params.m(max(i,j):data.params.n))*data.params.l(i)*data.params.l(j);

data.freq = 0.5;

%%
global p
p = @progressbar;
p()
%% Finite-time LQR
% terminal Gain
% -------------
% PT = 0.01*eye(3*(4+2*n));
% xf = reshape(PT,numel(PT),1);
% disp('Calculation Riccati Gains...') ;
% odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% [time, P] = ode15s(@odefun_riccati, [10 0], xf, odeopts, data) ;
% save('gainMatrix_sQnlink.mat','P','time');

load('gainMatrix_sQnlink.mat');
data.lqr.P = P;
data.lqr.time = time;


%% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% [trajd0] = get_ini_cond(data.params);

% Zero Initial Error- Configuration
% ---------------------------------
% [trajd0] = Flat2state(get_agrresive_traj(0),data.params);
[trajd0] = get_trajd(-1,data);

xL0 = trajd0.xL(data.params.n).x ;
vL0 = trajd0.xL(data.params.n).dx{1};

R0 = trajd0.R;
Omega0 = trajd0.Omega;

q0 = [];
dq0 = [];
omega0 = [];
for i = 1:data.params.n
    q0 = [q0,trajd0.q(i).q];
    dq0 = [dq0, trajd0.q(i).dq{1}];
    omega0 = [omega0, vec_cross(trajd0.q(i).q,trajd0.q(i).dq{1})];
end

% state - structure
% -----------------
% [xL; vL; R; Omega; qi; omegai] - for = 1:n
% setting up x0 (initial state)
% -----------------------------
x0 = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(omega0,numel(omega0),1)];

%% SIMULATION
% ==========
disp('Simulating...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% odeopts = [] ;
[t, x] = ode15s(@odefun_singleQuadFlex, [0 10], x0, odeopts, data) ;

% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(0.1*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_,f_,M_] =  odefun_singleQuadFlex(t(i),x(i,:)',data);
   xd(i,:) = xd_';
   psi_exL(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evL(i) = norm(x(i,4:6)-xd(i,4:6));
   psi_R(i) = abs(0.5*trace(eye(3,3) - reshape(xd_(7:15),3,3)'*reshape(x(i,7:15),3,3)));
   f(i,1)= f_;
   M(i,:)= M_';
end

%% PLOTS
% ======
figure;
subplot(1,3,1);
plot(t(ind),psi_exL(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Position Error');
legend('psi-exL');
subplot(1,3,2);
plot(t(ind),psi_evL(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Velocity Error');
legend('psi-evL');
subplot(1,3,3);
plot(t(ind),psi_R(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Rotation Error');
legend('psi-R');

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
anime_sQnlink(t(ind), x(ind,:),xd(ind,:), data.params,0);

end

%%
function[trajd0] = get_ini_cond(params)
traj.x = zeros(3,1);

for i = 1:(4+2*params.n)
    traj.dx{i} = zeros(3,1);
end
   
[trajd0]= Flat2state(traj,params);
th = 2/180*pi;
unitvector = [sin(th);0;cos(th)];
for i = 1:params.n
   trajd0.q(i).q = unitvector;
endtrajd0.R =  RPYtoRot_ZXY(178*pi/180,0,0);
end

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

O = zeros(3);
I = eye(3);


% fetching desired states
% -----------------------
% [trajd] = get_desired_traj(get_flat_traj(t),data.params);
% [trajd] = Flat2state(get_agrresive_traj(t),data.params);
[trajd] = get_trajd(t,data);
xLd = trajd.xL(data.params.n).x ;
vLd = trajd.xL(data.params.n).dx{1};
aLd = trajd.xL(data.params.n).dx{2};
Rd = trajd.R;
Omegad = trajd.Omega;

fd = trajd.f;
Md = trajd.M;

qd = [];
dqd = [];
omegad = [];
for i = 1:data.params.n
    qd = [qd,trajd.q(i).q];
    dqd = [dqd, trajd.q(i).dq{1}];
    omegad = [omegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
    
end
xQd = xLd - sum(repmat(l,3,1).*qd,2);
vQd = vLd - sum(repmat(l,3,1).*dqd,2);

% Extracing states
% ----------------
xL = x(1:3);
vL = x(4:6);
R = reshape(x(7:15),3,3);
Omega = x(16:18);
q = reshape(x(19:3*n+18),3,n);
omega = reshape(x(3*n+19:6*n+18),3,n);
dq = zeros(size(q));
for i = 1:n
   dq(:,i) = vec_cross(omega(:,i),q(:,i));
end
b3 = R(:,3);

xQ = xL - sum(repmat(l,3,1).*q,2);
vQ = vL - sum(repmat(l,3,1).*dq,2);

% =========================================================================
% calculating error state 's'
s = [(xQ-xQd); (vQ-vQd);];   %single link different
eta = 0.5*vee(Rd'*R - R'*Rd);
del_Om = Omega - (R'*Rd)*Omegad;

xi = []; del_om=[];
for i = 1:n
    xi = [xi;vec_cross(qd(:,i),q(:,i))];
    del_om = [del_om; omega(:,i)+vec_cross(q(:,i),vec_cross(q(:,i) ,omegad(:,i) ))];
end
s = [s; eta; del_Om; xi; del_om];

gain_K = get_linear_control(t,trajd,data);

du = gain_K*s;

f = fd + du(1);
M = Md +du(2:4);

% =========================================================================
    % Equations of Motion
    % -------------------
    dx = [];

    xL_dot = vL;
    q_dot = dq;

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
function [K] = get_linear_control(t,trajd,data)
%
[A,B] = get_linear_error_dynamics(trajd,data);
% gainQ = 0.01*eye(18);
% gainR = 0.05*eye(4);
Q1 = blkdiag(0.5*eye(6), eye(3), 0.75*eye(9));
Q2 = 0.2*eye(4);
Q2i = inv(Q2);

P_ = interp1(data.lqr.time,data.lqr.P,t);
m = sqrt(numel(P_));
P = reshape(P_,m,m);

K = -Q2i*B'*P;

end

%%
function [A,B] = get_linear_error_dynamics(trajd,data)
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

O = zeros(3);
I = eye(3);

xLd = trajd.xL(data.params.n).x ;
vLd = trajd.xL(data.params.n).dx{1};
aLd = trajd.xL(data.params.n).dx{2};
Rd = trajd.R;
Omegad = trajd.Omega;

fd = trajd.f;
Md = trajd.M;

qd = [];
dqd = [];
d2qd = [];
omegad = [];
domegad = [];
for i = 1:data.params.n
    qd = [qd,trajd.q(i).q];
    dqd = [dqd, trajd.q(i).dq{1}];
    d2qd = [d2qd,trajd.q(i).dq{2}];
    omegad = [omegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
    domegad = [domegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{2})];
end

wd = omegad;
dwd = domegad;

xQd = xLd - sum(repmat(l(1:n),3,1).*qd,2);
vQd = vLd - sum(repmat(l(1:n),3,1).*dqd,2);
aQd = aLd - sum(repmat(l(1:n),3,1).*d2qd,2);
% =========================================================================
% VARIATION-BASED LINEARIZATION: n LINKS
% (E)x = (F)s + (G)du => x = (E\F)s + (E\G)du

E = M00*I;
E1 = [];E2 = [];E3 = zeros(3*n);
for i = 1:n
    E1 = [E1, -M0i(i)*hat(qd(:,i))];
    E2 = [E2; Mi0(i)*hat(qd(:,i))];
    for j = 1:n
        if j==i
            E3(3*i-2:3*i,3*j-2:3*j) = Mij(i,j)*I;
        else
            E3(3*i-2:3*i,3*j-2:3*j) = -Mij(i,j)*hat(qd(:,i))*hat(qd(:,j));
        end
    end
end
E = [E, E1; E2, E3];

F = [O O -fd*Rd*hat(e3) O;];
F = [F;zeros(3*n,3*4)];
F1 = [];F2=[];F3 = zeros(3*n);F4 = zeros(3*n);

for i = 1:n
    F1 = [F1, M0i(i)*(hat(dwd(:,i))-norm(wd(:,i))^2)*hat(qd(:,i))];
    F2 = [F2, M0i(i)*2*qd(:,i)*wd(:,i)'];
    for j=1:n
        if i == j
            temp = Mi0(i)*hat(aQd) + sum(m(i:n))*g*l(i)*hat(e3);
            for k = 1:n
                if k~=i
                   temp = temp - Mij(i,k)*(hat(hat(qd(:,k))*dwd(:,k))+norm(wd(:,k))^2*hat(qd(:,k))); 
                end
            end
            F3(3*i-2:3*i,3*j-2:3*j) = temp*(-hat(qd(:,i)));
            F4(3*i-2:3*i,3*j-2:3*j) = O;
        else
            F3(3*i-2:3*i,3*j-2:3*j) = Mij(i,j)*hat(qd(:,i))*(hat(dwd(:,j))-norm(wd(:,j))^2)*hat(qd(:,j));
            F4(3*i-2:3*i,3*j-2:3*j) = 2*Mij(i,j)*hat(qd(:,i))*qd(:,j)*wd(:,j)';
        end
    end
end
    
F = [F, [F1, F2; F3, F4]];

G = [Rd*e3, O; zeros(3*n,4)];    

F_ = E\F;
G_ = E\G;

A33 = [];
A34 = [];
for i = 1:n
   A33 = blkdiag(A33,qd(:,i)*qd(:,i)'*hat(wd(:,i)));
   A34 = blkdiag(A34,I-qd(:,i)*qd(:,i)');
end

A66 = inv(J)*(hat(J*Omegad) - hat(Omegad)*J);
B62 = inv(J);

A = [O I zeros(3,3*(2+2*n));
        F_(1:3,:);
        O O -hat(Omegad) I zeros(3,3*(2*n));
        O O O           A66 zeros(3,3*(2*n))];
A = [A;[zeros(3*n,3*4), A33, A34]];
A = [A;F_(4:end,:)];



B = [zeros(3,1) O;
        G_(1:3,:);
        zeros(3,1) O;
        zeros(3,1) B62]; 
B = [B;
        zeros(3*n,4);
        G_(4:end,:)];


end

%%
function [dx] = odefun_riccati(t,x,data)
% local parameters/gains
n = data.params.n;
Q1 = blkdiag(0.5*eye(6),0.75*eye(6), eye(3*n), 0.75*eye(3*n));

Q2 = 0.2*eye(4);
Q2i = inv(Q2);

m = sqrt(numel(x));
P = reshape(x,m,m);

[trajd] = get_trajd(t,data);
[A,B] = get_linear_error_dynamics(trajd,data);

dP = - (Q1 - (P*B*Q2i*B'*P) + A'*P + P*A);

dx = reshape(dP,numel(dP),1);

if nargout <= 1
   fprintf('Sim time %0.4f seconds \n',t);
end
% global p
% p(t/25);\

end

function [trajd] = get_trajd(t,data)

[trajd] = Flat2state(get_agrresive_traj(t),data.params);
% [trajd] = Flat2state(testTraj(t,plannedTraj),data.params);
end




