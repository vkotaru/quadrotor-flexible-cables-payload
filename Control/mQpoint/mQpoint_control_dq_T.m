function mQpoint_control_dq_T()
%% sQ_nlink_control
% by vkotaru@andrew.cm.edu
% Date: July-18-2016
% Last Updated: July-18-2016
% -----------------------------------------------------------------------
% control 
% -----------------------------------------------------------------------
% quad_plan_traj() 
%
% Input Argument: 
% 				None 
% Return:
% 		None
% -----------------------------------------------------------------------

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
data.nQ = 3; % No. of Quadrotors
data.mL = 0.75; % Mass of the suspended Load
data.g = 9.81 ;
data.e1 = [1;0;0] ;
data.e2 = [0;1;0] ;
data.e3 = [0;0;1] ;
data.freq = 0.5;

% ---- Quad-1 ----------%
data.params(1).mQ = 0.85 ;
data.params(1).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(1).g = 9.81 ;
data.params(1).e1 = [1;0;0] ;
data.params(1).e2 = [0;1;0] ;
data.params(1).e3 = [0;0;1] ;

data.params(1).m = [0.1,0]; % 0.1*ones(1,5); % mass of each link
data.params(1).l = 0.5*ones(1,2); % length of each link 
data.params(1).n = 2; % No. of links suspended

% ---- Quad-2 ----------%
data.params(2).mQ = 0.85 ;
data.params(2).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(2).g = 9.81 ;
data.params(2).e1 = [1;0;0] ;
data.params(2).e2 = [0;1;0] ;
data.params(2).e3 = [0;0;1] ;

data.params(2).m = [0.1,0]; % 0.1*ones(1,5); % mass of each link
data.params(2).l = 0.5*ones(1,2); % length of each link 
data.params(2).n = 2; % No. of links suspended

% ---- Quad-3 ----------%
data.params(3).mQ = 0.85 ;
data.params(3).J = diag([0.557, 0.557, 1.05]*10e-2);
data.params(3).g = 9.81 ;
data.params(3).e1 = [1;0;0] ;
data.params(3).e2 = [0;1;0] ;
data.params(3).e3 = [0;0;1] ;

data.params(3).m = [0.1,0]; % 0.1*ones(1,5); % mass of each link
data.params(3).l = 0.5*ones(1,2); % length of each link 
data.params(3).n = 2; % No. of links suspendeddata.params.mQ = 0.85 ;

% data.params.M00 = data.params.mQ + sum(data.params.m(1:data.params.n));
% data.params.M0i = @(i) sum(data.params.m(i:data.params.n))*data.params.l(i); 
% data.params.Mi0 = data.params.M0i;
% data.params.Mij = @(i,j) sum(data.params.m(max(i,j):data.params.n))*data.params.l(i)*data.params.l(j);

data.freq = 0.5;
%% Finite-time LQR
% terminal Gain
% -------------
% PT = 0.01*eye(60);
% xf = reshape(PT,numel(PT),1);
% disp('Calculation Riccati Gains...') ;
% odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% [time, P] = ode15s(@odefun_riccati, [20 0], xf, odeopts, data) ;
% save('mQ_gainMatrix.mat','P','time');


load('mQ_gainMatrix.mat');
data.lqr.P = P;
data.lqr.time = time;


%% INTIALIZING - INTIAL CONDITIONS
% ================================
% Zero Position 
% -------------
% [trajd0] = get_ini_cond(data);

% Zero Initial Error- Configuration
% ---------------------------------
[trajd0] = get_trajd(-0.5,data);
% [A,B] = get_linearDynamics_mQpt(trajd0,data);
x0 = [];
data.ns = zeros(data.nQ,1);
for k = 1:data.nQ
    n = data.params(k).n-1;
    xL0 = trajd0(k).xL(n).x ;
    vL0 = trajd0(k).xL(n).dx{1};

    R0 = trajd0(k).R;
    Omega0 = trajd0(k).Omega;

    q0 = [];
    dq0 = [];
    omega0 = [];
    for i = 1:n
        q0 = [q0,trajd0(k).q(i).q];
        dq0 = [dq0, trajd0(k).q(i).dq{1}];
        omega0 = [omega0, vec_cross(trajd0(k).q(i).q,trajd0(k).q(i).dq{1})];
    end
    x0_ = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(omega0,numel(omega0),1)];
    data.ns(k) = length(x0_);
    x0 = [x0; x0_];
end

xLoad = trajd0(1).xL(data.params(1).n).x;
vLoad = trajd0(1).xL(data.params(1).n).dx{1};

% state :- ([xL; vL; R; Omega; qi; omegai])j - for = 1:n, j=1:nQ
x0 = [xLoad;vLoad; x0];

%% SIMULATION
% ========================================================================
disp('Simulating...') ;
odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9) ;
% odeopts = [] ;
[t, x] = ode15s(@odefun_mQPointFlex, [0 5], x0, odeopts, data) ;
% ========================================================================
% Computing Various Quantities
disp('Computing...') ;
ind = round(linspace(1, length(t), round(0.1*length(t)))) ;
% ind = 0:length(t);
for i = ind
   [~,xd_,] =  odefun_mQPointFlex(t(i),x(i,:)',data);
   xd(i,:) = xd_';

   % Load Error
   psi_exLoad(i) = norm(x(i,1:3)-xd(i,1:3));
   psi_evLoad(i) = norm(x(i,4:6)-xd(i,4:6));
%    
   s = mat2cell(x(i,7:end),1,data.ns');
   s_d = mat2cell(xd(i,7:end),1,data.ns');
   for j = 1:data.nQ
        psi_ex{j}(i) = norm(s{j}(1,1:3)-s_d{j}(1,1:3));
        psi_ev{j}(i) = norm(s{j}(1,4:6)-s_d{j}(1,4:6));
        psi_R{j}(i) = abs(0.5*trace(eye(3,3) - reshape(s_d{j}(1,7:15),3,3)'*reshape(s{j}(1,7:15),3,3)));
        lastairbender{j}(i) = norm(x(i,1:3)-s{j}(1,1:3));
   end
   
end

%% PLOTS
% ======
% figure;
% subplot(1,3,1);
% plot(t(ind),psi_exL(ind)); grid on;
% xlabel('time [s]');ylabel('error');title('Position Error');
% legend('psi-exL');
% subplot(1,3,2);
% plot(t(ind),psi_evL(ind)); grid on;
% xlabel('time [s]');ylabel('error');title('Velocity Error');
% legend('psi-evL');
% subplot(1,3,3);
% plot(t(ind),psi_R(ind)); grid on;
% xlabel('time [s]');ylabel('error');title('Rotation Error');
% legend('psi-R');
% 
figure;
subplot(1,2,1);
plot(t(ind),psi_exLoad(ind)); grid on;title('error x-Load');grid on;
subplot(1,2,2);
plot(t(ind),psi_exLoad(ind)); grid on;title('error v-Load');grid on;

for j = 1:data.nQ
    figure;
    subplot(2,2,1);
    plot(t(ind),psi_ex{j}(ind)); grid on;
    xlabel('time [s]');ylabel('error');title('Position Error');
    legend('psi-exL');
    subplot(2,2,2);
    plot(t(ind),psi_ev{j}(ind)); grid on;
    xlabel('time [s]');ylabel('error');title('Velocity Error');
    legend('psi-evL');
    subplot(2,2,3);
    plot(t(ind),psi_R{j}(ind)); grid on;
    xlabel('time [s]');ylabel('error');title('Rotation Error');
    legend('psi-R');
    subplot(2,2,4);
    plot(t(ind),lastairbender{j}(ind)); grid on;
    title('last link length');
end

keyboard;
anime_mQpt(t(ind), x(ind,:),xd(ind,:), data,0);

end

%%
function[dx, xd] = odefun_mQPointFlex(t,x,data)
% parameters
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1; 
e2 = data.e2;
e3 = data.e3; 
freq = data.freq;

% fetching desired states
[trajd] =  get_trajd(t,data);
xLoadd = trajd(1).xL(data.params(1).n).x;
vLoadd = trajd(1).xL(data.params(1).n).dx{1};

% extracting state
xLoad = x(1:3);
vLoad = x(4:6);
% rest of the states
x = x(7:end);
states = mat2cell(x',1,data.ns');
errorState = [];
for k = 1:nQ
    % desired state
    n = data.params(k).n-1;
    l = data.params(k).l;
    xLd = trajd(k).xL(n).x ;
    vLd = trajd(k).xL(n).dx{1};
    Rd = trajd(k).R;
    Omegad = trajd(k).Omega;
    qd = [];
    dqd = [];
    omegad = [];
    for i = 1:n
        qd = [qd,trajd(k).q(i).q];
        dqd = [dqd, trajd(k).q(i).dq{1}];
        omegad = [omegad, vec_cross(trajd(k).q(i).q,trajd(k).q(i).dq{1})];
    end
    xQd = trajd(k).xQ.x;vQd = trajd(k).xQ.dx{1};
    % actual state
    x = states{k}';
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
    
    % error calculation
    s = [(xQ-xQd); (vQ-vQd);];   %single link different
    eta = 0.5*vee(Rd'*R - R'*Rd);
    del_Om = Omega - (R'*Rd)*Omegad;

    xi = []; del_om=[];
    for i = 1:n
        xi = [xi;vec_cross(qd(:,i),q(:,i))];
        del_om = [del_om; omega(:,i)+vec_cross(q(:,i),vec_cross(q(:,i) ,omegad(:,i) ))];
    end
s = [eta; del_Om; (xQ-xQd); (vQ-vQd); xi; del_om]; 
errorState = [errorState;s]  ;
end
% ========================================================================
% CONTROL   
errorState = [errorState;(xLoad-xLoadd); (vLoad-vLoadd)];
% 
[gain_K] = get_linear_control(t,trajd,data);

du = gain_K*errorState;

du = reshape(du,4,nQ);
% ========================================================================
% ===== Equations of Motion =====

Right = []; Left_= cell(1,data.nQ);
L1 = []; L2 = [];
sizesofdx = [];
for i = 1:nQ
    [Left_{i}, Right_, L1_, L2_,info(i)] = get_quad_dynamics(states{i}', xLoad, vLoad, data, trajd(i).f+du(1,i),trajd(i).M+du(2:4,i) ,i);
    Right = [Right;Right_];
    L1 = [L1; L1_];
    L2 = [L2, L2_];
    sizesofdx = [sizesofdx, 3*(info(i).n+1)+1];
end
Left = blkdiag(Left_{:});
Left = [Left L1; L2 eye(3)];
Right = [Right; -g*e3];

d2x = Left\Right;
vLoad_dot  = d2x(end-2:end);

dstates = mat2cell(d2x(1:end-3)',1,sizesofdx);
dx = [];
for i = 1:nQ
    dx_ = dstates{i}';
    dxL = info(i).vL;
    q = info(i).q;
    dq = info(i).dq;
    l = info(i).l;
    
    dvQ = dx_(1:3);
    d2q = reshape(dx_(4:end-1),3,info(i).n);
    
    domega = zeros(size(d2q));
    dvL = dvQ;
    for j = 1:info(i).n
        dvL = dvL + l(j)*d2q(:,j);
        domega(:,j) = vec_cross(q(:,j),d2q(:,j));
    end
    
    dx_ = [dxL; dvL; reshape(info(i).dR,9,1); info(i).dOmega;...
        reshape(dq,numel(dq),1); reshape(domega,numel(domega),1)];
    dx = [dx; dx_];
end


% Computing xd
% ------------
xd = [];
for k = 1:data.nQ
    n = data.params(k).n-1;
    xLd = trajd(k).xL(n).x ;
    vLd = trajd(k).xL(n).dx{1};

    Rd = trajd(k).R;
    Omegad = trajd(k).Omega;

    qd = [];
    dqd = [];
    omegad = [];
    for i = 1:n
        qd = [qd,trajd(k).q(i).q];
        dqd = [dqd, trajd(k).q(i).dq{1}];
        omegad = [omegad, vec_cross(trajd(k).q(i).q,trajd(k).q(i).dq{1})];
    end
    xd_ = [xLd; vLd; reshape(Rd,9,1); Omegad; reshape(qd,numel(qd),1); reshape(omegad,numel(omegad),1)];
    xd = [xd; xd_];
end
xLoadd = trajd(1).xL(data.params(1).n).x;
vLoadd = trajd(1).xL(data.params(1).n).dx{1};
xd = [xLoadd;vLoadd; xd];

% Computing dx
%-------------
dx = [vLoad;
      vLoad_dot;dx];

if nargout <= 1
   fprintf('Sim time %0.4f seconds \n',t);
end
    
end


function[LeftMatrix,RightMatrix,L1,L2, state_info] = get_quad_dynamics(x, xLoad, vLoad,data, f, M, qNo)
% Extracing parameters
% --------------------
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1; 
e2 = data.e2;
e3 = data.e3; 

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

O = zeros(3);
I = eye(3);

% fetching desired states
% -----------------------
% f = trajd.f;
% M = trajd.M;
% 
% Td = -trajd.B(end).B;

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

xQ = xL - sum(repmat(l(1:n),3,1).*q,2);
vQ = vL - sum(repmat(l(1:n),3,1).*dq,2);

the_l = norm(xL-xLoad);
the_q = (xLoad-xL)/l(data.params(qNo).n);
the_dq = (vLoad-vL)/l(data.params(qNo).n);

% =========================================================================
    % Equations of Motion
    % -------------------
    dx = [];

%     xL_dot = vL;
    q_dot = zeros(size(q));
    for i = 1:n
       q_dot(:,i) = dq(:,i);%vec_cross(omega(:,i),q(:,i));
    end

%     lhsMat = [];
    lhsMat = M00*eye(3);
    tmp_l = zeros(3*n,3*n);
    tmp_l2 = [];
    tmp_l3 = [-I*the_q];
    tmp_b1 = -(xLoad-xL)';
    for ii = 1:n
        lhsMat = [lhsMat, M0i(ii)*eye(3)];
        tmp_l2 = [tmp_l2; -hat_map(q(:,ii))^2*Mi0(ii)];
        tmp_l3 = [tmp_l3; l(ii)*hat(q(:,ii))^2*the_q];
        tmp_b1 = [tmp_b1, -l(ii)*(xLoad-xL)'];
        for jj = 1:n
            if ii==jj
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = Mij(ii,jj)*eye(3);
            else
                tmp_l(3*ii-2:3*ii,3*jj-2:3*jj) = -Mij(ii,jj)*hat_map(q(:,ii))^2;
            end
        end
    end
    tmp_b1 = [tmp_b1, 0];
    lhsMat = [lhsMat;tmp_l2, tmp_l];
    lhsMat = [lhsMat, tmp_l3];
    lhsMat = [lhsMat; tmp_b1];
    
    L1 = [zeros(3*(n+1),3);(xLoad-xL)'];
    L2 = [zeros(3,3*(n+1)) I*the_q/mL];
    
%     rhsMat = [];
    rhsMat = f*R*e3 - M00*g*e3;
    for i = 1:n
        tmp_rhs = -norm(dq(:,i))^2*Mij(i,i)*q(:,i) + sum(m(i:n))*g*l(i)*hat_map(q(:,i))^2*e3;
        rhsMat = [rhsMat; tmp_rhs];
    end
    rhsMat = [rhsMat; -norm(vLoad-vL)^2];
    
    LeftMatrix = lhsMat; RightMatrix = rhsMat;
    
    % calculated information
    state_info.n = n;
    state_info.vL = vL;
    state_info.q = q;
    state_info.dq = dq;
    state_info.omega = omega;
    state_info.the_q = the_q;
    state_info.l = l(1:n);
    
    % Quadrotor Attitude
    state_info.dR = R*hat_map(Omega) ;
    state_info.dOmega = (J)\( -vec_cross(Omega, J*Omega) + M ) ;    
end

%% LINEARIZATION
function [A,B] = get_linearDynamics_mQpt(desired_trajectory,data)
% Extracing parameters
% --------------------
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1; 
e2 = data.e2;
e3 = data.e3; 
freq = data.freq;

Efinal = [];Es = [];Eb = [];
Ffinal = [];Fs = [];Fb = []; Fe = zeros(3,6);
Gfinal = [];

for qNo = 1:nQ
    trajd = desired_trajectory(qNo);
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

    Td = trajd.B(data.params(qNo).n).norm;
    the_qd = trajd.q(data.params(qNo).n).q;
    L = l(data.params(qNo).n);
    % =========================================================================
    % VARIATION-BASED LINEARIZATION: n LINKS
    % (E)x = (F)s + (G)du => x = (E\F)s + (E\G)du

    E = [M00*I -M0i(1)*hat(qd) -(xLd-x1d)/L;
            hat(qd)*Mi0(1) Mij(1,1)*I -l(1)*hat(qd)*(xLd-x1d)/L;
            -(xLd-x1d)' l(1)*(xLd-x1d)'*hat(qd) 0 ];
        Eside = [O;O;(xLd-x1d)'];
        Ebottom = [O O (xLd-x1d)/(mL*L)];

    F = [-fd*Rd*hat(e3) O -Td*I/L O (M0i(1)*(hat(dwd)-norm(wd)^2)+Td*l(1)/L)*hat(qd) (2*M0i(1)*qd*wd') ;
            O O -Td*l(1)*hat(qd)/L O (Td*l(1)^2*hat(qd)^2/L + Td*l(1)*hat(xLd-x1d)*hat(qd)/L - m(1)*g*l(1)*hat(e3)*hat(qd)) O ;
            zeros(1,3) zeros(1,3) (aLd-a1d)' 2*(vLd-v1d)' (-(aLd-a1d)'*l(1)*hat(qd)+l(1)*(xLd-x1d)'*(hat(dwd)-norm(wd)^2)*hat(qd)-2*l(1)*(vLd-v1d)'*hat(wd)*hat(qd)) -2*l(1)*((xLd-x1d)'*qd*wd'+(vLd-v1d)'*hat(qd)) ];
    Fside = [Td*I/L O; Td*l(1)*hat(qd)/L O; -(aLd-a1d)' -2*(vLd-v1d)'];
    Fbottom = [O O Td*I/(mL*L) O -Td*l(1)*hat(qd)/(mL*L) O];
    Fend = [-Td*I/(mL*L) O];
    G = [Rd(:,3) O;
            zeros(4,4)];

    Efinal = blkdiag(Efinal,E);
    Es = [Es;Eside];
    Eb = [Eb,Ebottom];
    
    Ffinal = blkdiag(Ffinal,F);
    Fs = [Fs;Fside];
    Fb = [Fb,Fbottom];
    Fe = Fe+Fend;
    Gfinal = blkdiag(Gfinal,G);
    
    %-------------------------------
    A66 = inv(J)*(hat(J*Omegad) - hat(Omegad)*J);
    B62 = inv(J);
    A33 = [];A34 = [];
    for i = 1:n
       A33 = blkdiag(A33,qd(:,i)*qd(:,i)'*hat(wd(:,i)));
       A34 = blkdiag(A34,I-qd(:,i)*qd(:,i)');
    end
    
    eso3{qNo} = [-hat(Omegad) I;O A66];
    es2{qNo} = [A33,A34];
    inso3{qNo}= [zeros(3,4);[zeros(3,1) B62]];
end

Esuper = [Efinal, Es; Eb, I];  
Fsuper = [Ffinal, Fs; Fb, Fe];
Gsuper = [Gfinal; zeros(3,4*nQ)];

F_ = inv(Esuper)*Fsuper;
G_ = inv(Esuper)*Gsuper;

mmm = size(F_,2);
A = [eso3{1}, zeros(6,mmm-6); 
        O O O I zeros(3,mmm-12);
        F_(1:3,:);
        O O O O es2{1} zeros(3,mmm-18);
        F_(4:6,:);
        zeros(6,18),eso3{2},zeros(6,mmm-24);
        zeros(3,24),O, I,zeros(3,mmm-30);
        F_(8:10,:);
        zeros(3,30),es2{2},zeros(3,mmm-36);
        F_(11:13,:);
        zeros(6,36),eso3{3},zeros(6,60-42);
        zeros(3,42),O,I,zeros(3,60-48);
        F_(15:17,:);
        zeros(3,48),es2{3},zeros(3,6);
        F_(18:20,:);
        zeros(3,54),O,I;
        F_(22:24,:)];

B = [inso3{1} zeros(6,8);
        zeros(3,12);
        G_(1:3,:);
        zeros(3,12);
        G_(4:6,:);
        zeros(6,4), inso3{2},zeros(6,4);
        zeros(3,12);
        G_(8:10,:);
        zeros(3,12);
        G_(11:13,:);
        zeros(6,8),inso3{3};
        zeros(3,12);
        G_(15:17,:);
        zeros(3,12);
        G_(18:20,:);
        zeros(3,12);
        G_(22:24,:);];
    
end

%%
function[trajd] = get_trajd(t,data)
% trajd = Flat2state_mQpoint(get_flat_traj_10link(t,data.freq),data);
trajd = Flat2state_mQpoint(get_agrresive_traj(t),data);
end

function [K] = get_linear_control(t,trajd,data)
%
[~,B] = get_linearDynamics_mQpt(trajd,data);
Q2 = 0.2*eye(12);
Q2i = inv(Q2);

P_ = interp1(data.lqr.time,data.lqr.P,t);
m = sqrt(numel(P_));
P = reshape(P_,m,m);

K = -Q2i*(B'*P);% + D'*C);

end

function [dx] = odefun_riccati(t,x,data)
% local parameters/gains
n = 1;
Q1 = blkdiag(0.5*eye(6),0.75*eye(6), eye(3*n), 0.75*eye(3*n),...
    0.5*eye(6),0.75*eye(6), eye(3*n), 0.75*eye(3*n),...
    0.5*eye(6),0.75*eye(6), eye(3*n), 0.75*eye(3*n),0.5*eye(6));
% Q1 = blkdiag(eye(6),eye(6), eye(3*n), eye(3*n),eye(6));

Q2 = 0.2*eye(12);
Q2i = inv(Q2);

m = sqrt(numel(x));
P = reshape(x,m,m);

[trajd] = get_trajd(t,data);

[A,B] = get_linearDynamics_mQpt(trajd,data);

dP = - (Q1 - (P*B*Q2i*B'*P) + A'*P + P*A);

dx = reshape(dP,numel(dP),1);

if nargout <= 1
   fprintf('tLQR %0.4f seconds \n',t);
end

end

function[trajd0] = get_ini_cond(data)
if nargin == 2
    traj.x = varargin{1};
else
    traj.x = zeros(3,1);
end

for i = 1:(4+4)
    traj.dx{i} = zeros(3,1);
end
   
[trajd0]= Flat2state_mQpoint(traj,data);
% th = 2/180*pi;
% unitvector = [sin(th);0;cos(th)];
% for i = 1:params.n
%    trajd0.q(i).q = unitvector;
% end
% trajd0.R =  RPYtoRot_ZXY(178*pi/180,0,0);

end