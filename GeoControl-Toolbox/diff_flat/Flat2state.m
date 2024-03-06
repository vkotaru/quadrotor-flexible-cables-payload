function[trajd]= Flat2state(traj,params)
% Function to generate Desired trajectory and states for Quadrotor with
% suspended flexible cable i.e., n links
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 18-May-2016
% =====================================================================

% Geometric Control Toolbox
% addpath('~/kotaru/Quad_FlexCable/GeoControl-Toolbox/');

% Temperary variables/[arameters
% load('params.mat');
% traj = get_flat_traj(0);

% =============
% Extracting Parameters
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
xL(n).x = traj.x;
xL(n).dx = traj.dx;

% =============
B = struct([]);
q = struct([]);
% Tensions in nth link
B(n).B = m(n)*(xL(n).dx{2} + g*e3);
for i = 1:(imax-2)
    B(n).dB{i} = m(n)*xL(n).dx{i+2};
end


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
    
    % desired tensions
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
trajd.xQ = xQ;

% Quad Attitude
trajd.R = R;
trajd.Omega = Omega;
trajd.dOmega = dOmega;
trajd.dR = dR;
trajd.d2R = d2R;
 
% inputs
trajd.f = f;
trajd.M = M;

% Chain Position
trajd.xL = xL;
trajd.q = q;
trajd.B = B;


end




