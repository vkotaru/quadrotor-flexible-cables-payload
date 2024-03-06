function sQnlink_Tpd()
%% no HELP
%% INITIALZING WORKSPACE
% =====================
% profile on
% Clear workspace
% clear; 
close all; 
% clc;

% Add Paths
% Geometric Control Toolbox
% addpath('../../GeoControl-Toolbox/');
% addpath('../');

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

data.freq = 0.5;
%% Finite-time LQR
% terminal Gain
% -------------
n = data.params.n-1;
PT = 0.01*eye(3*(6+2*n));
xf = reshape(PT,numel(PT),1);
disp('Calculation Riccati Gains...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
[time, P] = ode15s(@odefun_riccati, [20 0], xf, odeopts, data) ;
save('sQ_nlink_lqrgains.mat','P','time');

% load('sQ_nlink_lqrgains.mat');
data.lqr.P = P;
data.lqr.time = time;


%% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% [trajd0] = get_ini_cond(data.params);

% Zero Initial Error- Configuration
% ---------------------------------
% [trajd0] = Flat2state(get_flat_traj_10link(-0.1,1),data.params);
% [trajd0] = Flat2state(get_agrresive_traj(-0.1),data.params);
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

% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(0.1*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_,T_,Td_] =  odefun_singleQuadFlex_dom(t(i),x(i,:)',data);
   xd(i,:) = xd_';
   T(i,:) = T_';
   Td(i,:) = Td_';
   
   % Load Error
   psi_exLoad(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evLoad(i) = norm(x(i,4:6)-xd(i,4:6));
   
   % last link length
   lastairbender(i) = norm(x(i,1:3) - x(i,7:9));
   
   % xn Error
   x1 = x(i,7:end); x2 = xd(i,7:end);
   psi_exL(i) = norm(x1(1:3)-x2(1:3));
   psi_evL(i) = norm(x1(4:6)-x2(4:6));
   psi_R(i) = (0.5*trace(eye(3,3) - reshape(x2(7:15),3,3)'*reshape(x1(7:15),3,3)));
   
   % Tensions Error
   errT(i) = norm(T(i,:)-Td(i,:));
   
end

%% PLOTS
% ======
figure;
subplot(2,3,1);
plot(t(ind),psi_exL(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Position Error');
legend('psi-exL');
subplot(2,3,2);
plot(t(ind),psi_evL(ind)); grid on;
xlabel('time [s]');ylabel('error');title('Velocity Error');
legend('psi-evL');
subplot(2,3,3);
plot(t(ind),psi_R(ind)); grid on;
xlabel('time [s]');ylabel('error');title('\Psi_R');
legend('psi-R');

subplot(2,3,4);
plot(t(ind),psi_exLoad(ind)); grid on;title('error x-Load');
subplot(2,3,5);
plot(t(ind),psi_evLoad(ind)); grid on;title('error v-Load');

% figure;
% plot(t(ind),T(ind,1),'r',t(ind),T(ind,2),'g',t(ind),T(ind,3),'b');hold on;grid on;
% plot(t(ind),Td(ind,1),':r',t(ind),Td(ind,2),':g',t(ind),Td(ind,3),':b','LineWidth',1);hold on;grid on;

figure;subplot(1,2,1);
plot(t(ind),lastairbender(ind));grid on;title('the last airbender');
subplot(1,2,2);
plot(t(ind),errT(ind));grid on;title('T difference');

keyboard;
anime_sQ_2link(t(ind), x(ind,:),xd(ind,:), data.params,0);

end

%% ODE FUNCTION
function[dx, xd,T,Td] = odefun_singleQuadFlex_dom(t,x,data)
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
for i = 1:data.params.n
    qd = [qd,trajd.q(i).q];
    dqd = [dqd, trajd.q(i).dq{1}];
    omegad = [omegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
    
end
xQd = xLd - sum(repmat(l(1:n),3,1).*qd(:,1:n),2);
vQd = vLd - sum(repmat(l(1:n),3,1).*dqd(:,1:n),2);

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

xQ = xL - sum(repmat(l(1:n),3,1).*q(:,1:n),2);
vQ = vL - sum(repmat(l(1:n),3,1).*dq(:,1:n),2);

the_l = norm(xL-xLoad);
the_q = (xLoad-xL)/l(data.params.n);
the_dq = (vLoad-vL)/l(data.params.n);
% q = [q, the_q];
% dq = [dq, the_dq];
% omega = [omega,vec_cross(the_q,the_dq)];
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
s = [(xLoad-xLoadd);(vLoad-vLoadd);s];
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% epsilon_bar = 0.8 ;
% kp_xy = 0.3/epsilon_bar^2 ; kd_xy = 0.6/epsilon_bar ;
% k1 = diag([kp_xy kp_xy 2]) ; k2 = diag([kd_xy kd_xy 1.5]) ;
% Tension = +k1*(xLoad-xLoadd)+k2*(vLoad-vLoadd);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gain_K = get_linear_control(t,trajd,data,Td);

du = gain_K*s;

f = fd + du(1);
M = Md +du(2:4);

% f = fd ;
% M = Md ;
% =========================================================================
    % Equations of Motion
    % -------------------
    dx = [];

    xL_dot = vL;
    q_dot = dq(:,1:n);%vec_cross(omega(:,i),q(:,i));]

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


%% LINEARIZED SYSTEM
function [A,B] = get_linearDynamics(trajd,data,Td)
Td = [0;0;0];
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

M00 = mQ + sum(m(1:n));
M0i = @(i) sum(m(i:n))*l(i); 
Mi0 = M0i;
Mij = @(i,j) sum(m(max(i,j):n))*l(i)*l(j);

O = zeros(3);
I = eye(3);

xLd = trajd.xL(n).x ;
vLd = trajd.xL(n).dx{1};
aLd = trajd.xL(n).dx{2};
Rd = trajd.R;
Omegad = trajd.Omega;

fd = trajd.f;
Md = trajd.M;

qd = [];
dqd = [];
d2qd = [];
omegad = [];
domegad = [];
for i = 1:n
    qd = [qd,trajd.q(i).q];
    dqd = [dqd, trajd.q(i).dq{1}];
    d2qd = [d2qd,trajd.q(i).dq{2}];
    omegad = [omegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
    domegad = [domegad, vec_cross(trajd.q(i).q,trajd.q(i).dq{2})];
end

wd = omegad;
dwd = domegad;

xQd = xLd - sum(repmat(l(1:n),3,1).*qd(:,1:n),2);
vQd = vLd - sum(repmat(l(1:n),3,1).*dqd(:,1:n),2);
aQd = aLd - sum(repmat(l(1:n),3,1).*d2qd(:,1:n),2);

mL = m(data.params.n);
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
            temp = Mi0(i)*hat(aQd) + sum(m(i:n))*g*l(i)*hat(e3)+ l(i)*hat(Td)*hat(qd(:,i));
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

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
B = [zeros(6,4);B];
    
epsilon_bar = 0.8 ;
kp_xy = 0.3/epsilon_bar^2 ; kd_xy = 0.6/epsilon_bar ;
k1 = diag([kp_xy kp_xy 2]) ; k2 = diag([kd_xy kd_xy 1.5]) ;

s = length(A);
tmpA = [O O; k1/mL k2/mL; zeros(3*(2+n),6)];
for i = 1:n
   tmpA = [tmpA; l(i)*hat(qd(:,i))*k1 l(i)*hat(qd(:,i))*k2]; 
end

A = [[O I;-k1/mL -k2/mL],zeros(6,s); tmpA, A];

end

%% NECESSARY FUNCTIONS
function [dx] = odefun_riccati(t,x,data)
% local parameters/gains
n = data.params.n-1;
% Q1 = blkdiag(0.5*eye(6),0.75*eye(6), eye(3*n), 0.75*eye(3*n),0.5*eye(6));
Q1 = blkdiag(0.5*eye(6*2),0.75*eye(6), eye(3*n), 0.75*eye(3*n));
Q2 = 0.2*eye(4);
Q2i = inv(Q2);

m = sqrt(numel(x));
P = reshape(x,m,m);

[trajd] = get_trajd(t,data);
[A,B] =get_linearDynamics(trajd,data, -trajd.B(end).B);

dP = - (Q1 - (P*B*Q2i*B'*P) + A'*P + P*A);

dx = reshape(dP,numel(dP),1);

if nargout <= 1
   fprintf('tLQR %0.4f seconds \n',t);
end

end

function [K] = get_linear_control(t,trajd,data,Td)
%
[A,B] = get_linearDynamics(trajd,data,Td);
Q2 = 0.2*eye(4);
Q2i = inv(Q2);

P_ = interp1(data.lqr.time,data.lqr.P,t);
m = sqrt(numel(P_));
P = reshape(P_,m,m);

K = -Q2i*B'*P;

end

function[trajd0] = get_ini_cond(params)
traj.x = zeros(3,1);

for i = 1:(4+2*params.n)
    traj.dx{i} = zeros(3,1);
end
   
[trajd0]= Flat2state(traj,params);
% th = 2/180*pi;
% unitvector = [sin(th);0;cos(th)];
% for i = 1:params.n
%    trajd0.q(i).q = unitvector;
% endtrajd0.R =  RPYtoRot_ZXY(90*pi/180,90*pi/180,0);
% end

end

%%
function[trajd] = get_trajd(t,data)
% [trajd] = Flat2state(get_flat_traj_10link(t,1),data.params);
[trajd] = Flat2state(get_agrresive_traj(t),data.params);
end