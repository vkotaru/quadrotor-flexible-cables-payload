function sQ_2link_control_tension()
%% INITIALZING WORKSPACE
% =====================
% profile on
% Clear workspace
% clear; 
close all; 
% clc;

% Add Paths
% Geometric Control Toolbox
addpath('../../GeoControl-Toolbox/');
addpath('../');

%% INITIALZING PARAMETERS
% ======================
% System constants and parameters
data.params.mQ = 0.85 ;
data.params.J = diag([0.557, 0.557, 1.05]*10e-2);
data.params.g = 9.81 ;
data.params.e1 = [1;0;0] ;
data.params.e2 = [0;1;0] ;
data.params.e3 = [0;0;1] ;

data.params.m = [0.5, 0.5]; % mass of each link
data.params.l = 0.5*ones(1,2); % length of each link 
data.params.n = 2; % 

%% Finite-time LQR
% terminal Gain
% -------------
n = data.params.n-1;
PT = 0.01*eye(3*(6+2*n));
xf = reshape(PT,numel(PT),1);
disp('Calculation Riccati Gains...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
[time, P] = ode15s(@odefun_riccati, [20 0], xf, odeopts, data) ;
% [time, P] = ODEeuler(@odefun_riccati, [20 0], xf, 0.001, data) ;
% save('gainMatrix.mat','P','time');

% load('gainMatrix.mat');
data.lqr.P = P;
data.lqr.time = time;


%% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% [trajd0] = get_ini_cond(data.params);

% Zero Initial Error- Configuration
% ---------------------------------
% [trajd0] = Flat2state(get_flat_traj_10link(-0.2,1),data.params);
% [trajd0] = Flat2state(get_agrresive_traj(-0.2),data.params);
[trajd0] = get_ini_cond(data.params);
n = data.params.n - 1;

xL0 = trajd0.xL(n).x ;
vL0 = trajd0.xL(n).dx{1};

R0 = trajd0.R;
Omega0 = trajd0.Omega;

q0 = [];
dq0 = [];
omega0 = [];
for i = 1:n
    q0 = [q0,trajd0.q(i).q];
    dq0 = [dq0, trajd0.q(i).dq{1}];
    omega0 = [omega0, vec_cross(trajd0.q(i).q,trajd0.q(i).dq{1})];
end

% state :- [xL; vL; R; Omega; qi; omegai] - for = 1:n
x0 = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(omega0,numel(omega0),1)];

xLoad0 = trajd0.xL(data.params.n).x ;
vLoad0 = trajd0.xL(data.params.n).dx{1};
x0 = [xLoad0; vLoad0; x0];
% 
% [A,B] = get_linear_error_dynamics_T(trajd0,data);

%% SIMULATION
% ==========
disp('Simulating...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% odeopts = [] ;
[t, x] = ode15s(@odefun_singleQuadFlex_dom, [0 20], x0, odeopts, data) ;
% [t, x] = ODEeuler(@odefun_singleQuadFlex_dom, [0 2], x0, 0.001, data) ;

% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(.1*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_,T_,Td_,dT_,f_,M_] =  odefun_singleQuadFlex_dom(t(i),x(i,:)',data);
   xd(i,:) = xd_';
   T(i,:) = T_';
   Td(i,:) = Td_';
   dT(i,:) = dT_';
   %forces
   f(i) = f_';
   M(i,:) = M_';
   % Load Error
   psi_exLoad(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evLoad(i) = norm(x(i,4:6)-xd(i,4:6));
   
   % last link length
   lastairbender(i) = norm(x(i,1:3) - x(i,7:9));
   
   % xn Error
   x1 = x(i,7:end); x2 = xd(i,7:end);
   psi_exL(i) = norm(x1(1:3)-x2(1:3));
   psi_evL(i) = norm(x1(4:6)-x2(4:6));
   psi_R(i) = abs(0.5*trace(eye(3,3) - reshape(x2(7:15),3,3)'*reshape(x1(7:15),3,3)));
   
   % Tensions Error
   errT(i) = norm(T(i,:)-Td(i,:));
   
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
subplot(1,2,1);
plot(t(ind),psi_exLoad(ind)); grid on;title('error x-Load');
subplot(1,2,2);
plot(t(ind),psi_evLoad(ind)); grid on;title('error v-Load');

figure;
plot(t(ind),T(ind,1),'r',t(ind),T(ind,2),'g',t(ind),T(ind,3),'b');hold on;grid on;
plot(t(ind),Td(ind,1),':r',t(ind),Td(ind,2),':g',t(ind),Td(ind,3),':b','LineWidth',1);hold on;grid on;
title('Tensions');
figure;
plot(t(ind),lastairbender(ind));grid on;title('the last airbender');

figure;
plot(t(ind),errT(ind),':r',t(ind),dT(ind),'b','LineWidth',1);grid on;title('difference in T');
legend('errorT','deltaT');

figure;
subplot(1,2,1);
plot(t(ind),f(ind),'b','LineWidth',1);grid on;title('thrust');
subplot(1,2,2);
plot(t(ind),M(ind,1),'r',t(ind),M(ind,2),'g',t(ind),M(ind,3),'b');grid on;
title('Moment');

keyboard;
anime_sQ_2link(t(ind), x(ind,:),xd(ind,:), data.params,0);

end


%%
function[dx, xd,T,Td,deltaT,f,M] = odefun_singleQuadFlex_dom(t,x,data)
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
n = data.params.n-1; % No. of links suspended

% M00 = data.params.M00;
% M0i = data.params.M0i; 
% Mi0 = data.params.Mi0;
% Mij = data.params.Mij;

M00 = mQ + sum(m(1:n));
M0i = @(i) sum(m(i:n))*l(i); 
Mi0 = M0i;
Mij = @(i,j) sum(m(max(i,j):n))*l(i)*l(j);

O = zeros(3);
I = eye(3);

mL = data.params.m(data.params.n);

% fetching desired states
% -----------------------
[trajd] = get_trajd(t,data);
xLd = trajd.xL(n).x ;
vLd = trajd.xL(n).dx{1};
aLd = trajd.xL(n).dx{2};
Rd = trajd.R;
Omegad = trajd.Omega;

fd = trajd.f;
Md = trajd.M;

qd = [];
dqd = [];
omegad = [];
for i = 1:n
    qd = [qd,trajd.q(i).q];
    dqd = [dqd, trajd.q(i).dq{1}];
    omegad = [omegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
    
end
xQd = trajd.xQ.x;
vQd = trajd.xQ.dx{1};

xLoadd = trajd.xL(data.params.n).x;
vLoadd = trajd.xL(data.params.n).dx{1};

Td = -trajd.B(end).B;

% Extracing states
% ----------------
xLoad = x(1:3);
vLoad = x(4:6);

x = x(7:end);
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

xQ = xL - sum(repmat(l(1:n),3,1).*q,2);
vQ = vL - sum(repmat(l(1:n),3,1).*dq,2);

the_l = norm(xL-xLoad);
the_q = (xLoad-xL)/l(data.params.n);
the_dq = (vLoad-vL)/l(data.params.n);
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
s = [eta; del_Om; (xQ-xQd); (vQ-vQd); xi; del_om; (xLoad-xLoadd); (vLoad-vLoadd)];
% 
[gain_K,C,D] = get_linear_control(t,trajd,data);

du = gain_K*s;
deltaT = C*s + D*du;


f = fd + du(1);
M = Md +du(2:4);
% 
% f = fd ;
% M = Md ;
% =========================================================================
    % Equations of Motion
    % -------------------
    dx = [];

    xL_dot = vL;
    q_dot = dq;%vec_cross(omega(:,i),q(:,i));]

    %     lhsMat = [];
    lhsMat = M00*eye(3);
    tmp_l = zeros(3*n,3*n);
    tmp_l2 = [];
    tmp_l3 = [-I*the_q O];
    tmp_b1 = -(xLoad-xL)';
    for ii = 1:n
        lhsMat = [lhsMat, -M0i(ii)*hat_map(q(:,ii))];
        tmp_l2 = [tmp_l2; hat_map(q(:,ii))*Mi0(ii)];
        tmp_l3 = [tmp_l3; -l(ii)*vec_cross(q(:,ii),the_q) O];
        tmp_b1 = [tmp_b1, l(ii)*(xLoad-xL)'*hat(q(:,ii))];
        for jj = 1:n
            if ii==jj
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = Mij(ii,jj)*eye(3);
            else
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = -Mij(ii,jj)*hat_map(q(:,ii))*hat_map(q(:,jj));
            end
        end
    end
    tmp_b = [zeros(3,3*(n+1)) I*the_q/mL I];
    tmp_b1 = [tmp_b1, 0 (xLoad-xL)'];
    lhsMat = [lhsMat;tmp_l2, tmp_l];
    lhsMat = [lhsMat, tmp_l3];
    lhsMat = [lhsMat; tmp_b; tmp_b1];
    
%     rhsMat = [];
    rhsMat = f*R*e3 - M00*g*e3;
    for i = 1:n
       rhsMat = rhsMat + M0i(i)*norm(omega(:,i))^2*q(:,i);  
    end
    temp_r = zeros(3,1);
    for i = 1:n
        tmp_rhs =  - sum(m(i:n))*g*l(i)*hat_map(q(:,i))*e3;
        temp_r = temp_r + l(i)*norm(omega(:,i))^2*q(:,i);
        for j = 1:n
            if j~=i
               tmp_rhs = tmp_rhs + Mij(i,j)*norm(omega(:,j))^2*hat_map(q(:,i))*q(:,j) ;
            end
        end
        rhsMat = [rhsMat; tmp_rhs];
    end
    rhsMat = [rhsMat; -g*e3; -(xLoad-xL)'*temp_r-norm(vLoad-vL)^2];
% ------------------------------------------------------------------------
    
    d2x = lhsMat\rhsMat;
    
    vLoad_dot = d2x(end-2:end);
    F = d2x(3*(n+1)+1);
    T = F*the_q;
    
    omega_dot = reshape(d2x(4:end-4),3,n);
    
    vL_dot = d2x(1:3,1);
    for i = 1:n
        vL_dot = vL_dot + l(i)*(vec_cross(omega_dot(:,i),q(:,i))+ vec_cross(omega(:,i),q_dot(:,i)));
    end
        
    % Quadrotor Attitude
    R_dot = R*hat_map(Omega) ;
    Omega_dot = (J)\( -vec_cross(Omega, J*Omega) + M ) ;
    
% Computing xd
% ------------
xd = [trajd.xL(data.params.n).x;trajd.xL(data.params.n).dx{1}];
xd = [xd;trajd.xL(n).x;trajd.xL(n).dx{1}];
xd = [xd;reshape(trajd.R, 9,1);trajd.Omega];
for i = 1:n
    xd = [xd;trajd.q(i).q];
end
for i = 1:n
    xd = [xd;vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
end
% Computing dx
%-------------
dx = [vLoad;
      vLoad_dot;
      xL_dot;
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
function [A,B, C,D] = get_linear_error_dynamics_T(trajd,data)
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
n = data.params.n-1; % No. of links suspended
mL = m(data.params.n);

M00 = mQ + sum(m(1:n));
M0i = @(i) sum(m(i:n))*l(i); 
Mi0 = M0i;
Mij = @(i,j) sum(m(max(i,j):n))*l(i)*l(j);
O = zeros(3);
I = eye(3);
Rd = trajd.R;
Omegad = trajd.Omega;

fd = trajd.f;
Md = trajd.M;

qd = [];
dqd = [];
d2qd = [];
wd = [];
dwd = [];
for i = 1:n
    qd = [qd,trajd.q(i).q];
    dqd = [dqd, trajd.q(i).dq{1}];
    d2qd = [d2qd,trajd.q(i).dq{2}];
    wd = [wd, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
    dwd = [dwd, vec_cross(trajd.q(i).q,trajd.q(i).dq{2})];
end
xQd = trajd.xQ.x ;
vQd = trajd.xQ.dx{1};
aQd = trajd.xQ.dx{2};

x1d = trajd.xL(1).x ;
v1d = trajd.xL(1).dx{1};
a1d = trajd.xL(1).dx{2};

xLd = trajd.xL(2).x ;
vLd = trajd.xL(2).dx{1};
aLd = trajd.xL(2).dx{2};

Td = trajd.B(data.params.n).norm;
the_qd = trajd.q(data.params.n).q;
L = l(data.params.n);
% =========================================================================
% VARIATION-BASED LINEARIZATION: n LINKS
% (E)x = (F)s + (G)du => x = (E\F)s + (E\G)du

E = [M00*I -M0i(1)*hat(qd) -(xLd-x1d)/L O;
        hat(qd)*Mi0(1) Mij(1,1)*I -l(1)*hat(qd)*(xLd-x1d)/L O;
        O O (xLd-x1d)/(mL*L) I;
        -(xLd-x1d)' l(1)*(xLd-x1d)'*hat(qd) 0 (xLd-x1d)'];

F = [-fd*Rd*hat(e3) O -Td*I/L O (M0i(1)*(hat(dwd)-norm(wd)^2)+Td*l(1)/L)*hat(qd) (2*M0i(1)*qd*wd') Td*I/L O;
        O O -Td*l(1)*hat(qd)/L O (Td*l(1)^2*hat(qd)^2/L + Td*l(1)*hat(xLd-x1d)*hat(qd)/L - m(1)*g*l(1)*hat(e3)*hat(qd)) O Td*l(1)*hat(qd)/L O;
        O O Td*I/(mL*L) O -Td*l(1)*hat(qd)/(mL*L) O -Td*I/(mL*L) O;
        zeros(1,3) zeros(1,3) (aLd-a1d)' 2*(vLd-v1d)' (-(aLd-a1d)'*l(1)*hat(qd)+l(1)*(xLd-x1d)'*(hat(dwd)-norm(wd)^2)*hat(qd)-2*l(1)*(vLd-v1d)'*hat(wd)*hat(qd)) -2*l(1)*((xLd-x1d)'*qd*wd'+(vLd-v1d)'*hat(qd)) -(aLd-a1d)' -2*(vLd-v1d)'];

G = [Rd(:,3) O;
        zeros(7,4)];
    
F_ = inv(E)*F;
G_ = inv(E)*G;

A66 = inv(J)*(hat(J*Omegad) - hat(Omegad)*J);
B62 = inv(J);A33 = [];A34 = [];
for i = 1:n
   A33 = blkdiag(A33,qd(:,i)*qd(:,i)'*hat(wd(:,i)));
   A34 = blkdiag(A34,I-qd(:,i)*qd(:,i)');
end

A = [-hat(Omegad) I zeros(3,18);
        O A66 zeros(3,18);
        O O O I O O O O;
        F_(1:3,:);
        O O O O A33 A34 O O;
        F_(4:6,:);
        zeros(3,18), O I;
        F_(8:10,:)];

B = [zeros(3,4); 
        [zeros(3,1) B62];
      zeros(3,4);
      G_(1:3,:);
      zeros(3,4);
      G_(4:6,:);
      zeros(3,4);
      G_(8:10,:);];
  
C = F_(7,:);
D = G_(7,:);
    
end


%%
function [dx] = odefun_riccati(t,x,data)
% local parameters/gains
n = data.params.n-1;
Q1 = blkdiag(0.5*eye(6),0.75*eye(6), eye(3*n), 0.75*eye(3*n),0.5*eye(6));
% Q1 = blkdiag(eye(6),eye(6), eye(3*n), eye(3*n),eye(6));

Q2 = 0.2*eye(4);
Q2i = inv(Q2);

m = sqrt(numel(x));
P = reshape(x,m,m);

[trajd] = get_trajd(t,data);

[A,B,~,~] = get_linear_error_dynamics_T(trajd,data);

dP = - (Q1 - (P*B*Q2i*B'*P) + A'*P + P*A);

dx = reshape(dP,numel(dP),1);

if nargout <= 1
   fprintf('tLQR %0.4f seconds \n',t);
end

end

%%
function[trajd0] = get_ini_cond(params,varargin)
if nargin == 2
    traj.x = varargin{1};
else
    traj.x = zeros(3,1);
end

for i = 1:(4+2*params.n)
    traj.dx{i} = zeros(3,1);
end
   
[trajd0]= Flat2state(traj,params);
% th = 2/180*pi;
% unitvector = [sin(th);0;cos(th)];
% for i = 1:params.n
%    trajd0.q(i).q = unitvector;
% end
% trajd0.R =  RPYtoRot_ZXY(178*pi/180,0,0);

end

%%
function [K,C,D] = get_linear_control(t,trajd,data)
%
[~,B,C,D] = get_linear_error_dynamics_T(trajd,data);
Q2 = 0.2*eye(4);
Q2i = inv(Q2);

P_ = interp1(data.lqr.time,data.lqr.P,t);
m = sqrt(numel(P_));
P = reshape(P_,m,m);

K = -Q2i*(B'*P);% + D'*C);

end


%% getting desired traj
function[trajd] = get_trajd(t,data)
% [trajd] = Flat2state(get_flat_traj_10link(t,1),data.params);
[trajd] = Flat2state(get_agrresive_traj(t),data.params);
% [trajd] = get_ini_cond(data.params,[0.5;0;0]);
end


