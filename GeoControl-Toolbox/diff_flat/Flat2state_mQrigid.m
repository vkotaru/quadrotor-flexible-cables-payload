function[trajd]= get_desired_traj_multiQuad_rigidload(traj,data)
% Multiple Quadrotor system - with suspended Point Mass load via n-link
% flexible cables
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 26-May-2016
% =====================================================================

% Geometric Control Toolbox
addpath('../GeoControl-Toolbox/');

% Temperary variables/Parameters
% load('params.mat');
% traj = get_flat_traj(0);
% =============

% Extracting Parameters
% =====================
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1;
e2 = data.e2;
e3 = data.e3;
r = data.r;
JL = data.JL;

% Load Orientations
traj.RL = RotX(5*pi/180);%eye(3);
traj.dRL = zeros(3,3);
traj.d2RL = zeros(3,3);
traj.OmegaL = vee_map(traj.RL'*traj.dRL) ;
traj.dOmegaL = vee_map( traj.dRL'*traj.dRL + traj.RL'*traj.d2RL) ;


% Tensions - Flat-Outputs
LAMBDA = [(1/4)*data.mL*data.g*vec_dot(RotY(pi/6)*data.e3,data.e1); 
                (1/4)*data.mL*data.g*vec_dot(RotY(pi/6)*data.e3,data.e2);
                (1/4)*data.mL*data.g*vec_dot(RotX(pi/6)*data.e3,data.e1); 
                (1/4)*data.mL*data.g*RotY(-pi/6)*data.e3];

xLcom.x = traj.x;
xLcom.dx = traj.dx;

W = [traj.RL'*mL*(xLcom.dx{2}+g*e3); JL*traj.dOmegaL + vec_cross(traj.OmegaL,JL*traj.OmegaL)];
            
imax = 4 + 2*max([data.params(:).n]);
traj.T = data.Phi_pinv*W + data.N*LAMBDA;
for i = 1:(imax-2)
   traj.dT{i} = data.Phi_pinv*[traj.RL'*mL*(xLcom.dx{i+2});zeros(3,1)]; 
end
            

for i = 1:nQ
    
    params = data.params(i);
    trajd(i) = get_singleQuad_traj(traj,data,params,i);
    
end

for i = 1:nQ
   trajd(i).xLcom = xLcom; 
end

end


function[trajectory] = get_singleQuad_traj(traj,data,params,quadNum)

% Extracting Parameters
% =====================
nQ = data.nQ; % No. of Quadrotors
mL = data.mL; % Mass of the suspended Load
g = data.g;
e1 = data.e1;
e2 = data.e2;
e3 = data.e3;
r = data.r;
JL = data.JL;

n = params.n; % No. of links
m = params.m; % Mass of links
l = params.l; % Length of links
mQ = params.mQ; % Quad-rotor mass
g = params.g ; 
e1 = params.e1 ;
e2 = params.e2 ;
e3 = params.e3 ;
J = params.J;

imax = 4+(2*n);
% =============

% Extracting Flat-Outputs
xLcom.x = traj.x;
xLcom.dx = traj.dx;
RL = traj.RL;
OmegaL = traj.OmegaL;

xL(n).x = xLcom.x + RL*r(:,quadNum);
xL(n).dx = xLcom.dx;

B(n).B = RL*traj.T(3*quadNum-2:3*quadNum,1);
for i = 1:(imax-2)
    B(n).dB{i} = RL*traj.dT{i}(3*quadNum-2:3*quadNum,1);
end
    
% =============
% B = struct([]);
q = struct([]);
% % Tensions in nth link
% B(n).B = m(n)*(xL(n).dx{2} + g*e3);
% for i = 1:(imax-2)
%     B(n).dB{i} = m(n)*xL(n).dx{i+2};
% end


normBn = norm(B(n).B);
B(n).norm = normBn;
% nth Unit vector and its higher derivatives
q(n).q = -B(n).B/normBn;
q(n).dq = get_q_derivaties_10link(B(n).B,normBn,B(n).dB{:});


% ===================================
% Expressions for all remaining links
for i = n-1:-1:1
    
    % desired positions
    xL(i).x = xL(i+1).x - l(i+1)*q(i+1).q;
    for j = 1:length(q(i+1).dq)
        xL(i).dx{j} = xL(i+1).dx{j} - l(i+1)*q(i+1).dq{j};
    end
    
    % desried tensions
    B(i).B = m(i)*(xL(i).dx{2}+g*e3) + B(i+1).B;
    B(i).norm = norm(B(i).B);
    for j = 1:(length(q(i+1).dq)-2)
        B(i).dB{j} = m(i)*xL(i).dx{j+2} + B(i+1).dB{j};
    end

    
    % desired attitudes
    q(i).q = -B(i).B/B(i).norm;
    q(i).dq = get_q_derivaties_10link(B(i).B,B(i).norm,B(i).dB{:});
    
end

% Quad-desired positions
% ======================
xQ.x = xL(1).x - l(1)*q(1).q;
for j = 1:length(q(1).dq)
    xQ.dx{j} = xL(1).dx{j} - l(1)*q(1).dq{j};
end
    
% Quad-desired Orientation
% ========================
% Extract from get_nom_traj.m
% by Dr.Koushil Sreenath

b1d = e1 ;
db1d = zeros(3,1) ;
d2b1d = zeros(3,1) ;

fb3 = mQ*(xQ.dx{2}+g*e3) + B(1).B;
norm_fb3 = norm(fb3) ;
f = norm_fb3 ;
b3 = fb3 / norm_fb3 ;
b3_b1d = vec_cross(b3, b1d) ;
norm_b3_b1d = norm(b3_b1d) ;
b1 = - vec_cross(b3, b3_b1d) / norm_b3_b1d ;
b2 = vec_cross(b3, b1) ;
R = [b1 b2 b3] ;

dfb3 = mQ*(xQ.dx{3}) + B(1).dB{1} ;
dnorm_fb3 = vec_dot(fb3, dfb3) / norm_fb3 ;
db3 = (dfb3*norm_fb3 - fb3*dnorm_fb3) / norm_fb3^2 ;
db3_b1d = vec_cross(db3, b1d) + vec_cross(b3, db1d) ;
dnorm_b3_b1d = vec_dot(b3_b1d, db3_b1d) / norm_b3_b1d ;
db1 = (-vec_cross(db3,b3_b1d)-vec_cross(b3,db3_b1d) - b1*dnorm_b3_b1d) / norm_b3_b1d ;
db2 = vec_cross(db3, b1) + vec_cross(b3, db1) ;
dR = [db1 db2 db3] ;
Omega = vee_map(R'*dR) ;

d2fb3 = mQ*(xQ.dx{4}) + B(1).dB{2};
d2norm_fb3 = (vec_dot(dfb3, dfb3)+vec_dot(fb3, d2fb3) - dnorm_fb3*dnorm_fb3) / norm_fb3 ;
d2b3 = ( (d2fb3*norm_fb3+dfb3*dnorm_fb3 - dfb3*dnorm_fb3-fb3*d2norm_fb3)*norm_fb3^2 - db3*norm_fb3^2*2*norm_fb3*dnorm_fb3 ) / norm_fb3^4 ;
d2b3_b1d = vec_cross(d2b3, b1d)+vec_cross(db3, db1d) + vec_cross(db3, db1d)+vec_cross(b3, d2b1d) ;
d2norm_b3_b1d = ( (vec_dot(db3_b1d,db3_b1d)+vec_dot(b3_b1d,d2b3_b1d))*norm_b3_b1d - vec_dot(b3_b1d, db3_b1d)*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
d2b1 = ( (-vec_cross(d2b3,b3_b1d)-vec_cross(db3,db3_b1d) - vec_cross(db3,db3_b1d)-vec_cross(b3,d2b3_b1d) - db1*dnorm_b3_b1d-b1*d2norm_b3_b1d )*norm_b3_b1d - db1*norm_b3_b1d*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
d2b2 = vec_cross(d2b3, b1)+vec_cross(db3, db1) + vec_cross(db3, db1)+vec_cross(b3, d2b1) ;
d2R = [d2b1 d2b2 d2b3] ;
dOmega = vee_map( dR'*dR + R'*d2R ) ; %vee_map( dR'*dR + R'*d2R, true ) ;

M = J*dOmega + vec_cross(Omega, J*Omega) ;

% FINAL DESIRED TRAJECTORY
% ========================
% Quad position
trajectory.xQ = xQ;

% Quad Attitude
trajectory.R = R;
trajectory.Omega = Omega;
trajectory.dOmega = dOmega;
trajectory.dR = dR;
trajectory.d2R = d2R;
 
% inputs
trajectory.f = f;
trajectory.M = M;

% Chain Position
trajectory.xL = xL;
trajectory.q = q;
trajectory.B = B;

% Load Orientation
trajectory.RL = traj.RL;
trajectory.OmegaL = traj.OmegaL;

end

