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

n = data.params.n;
data.params.M00 = data.params.mQ + sum(data.params.m(1:data.params.n));
data.params.M0i = @(i) sum(data.params.m(i:data.params.n))*data.params.l(i); 
data.params.Mi0 = data.params.M0i;
data.params.Mij = @(i,j) sum(data.params.m(max(i,j):data.params.n))*data.params.l(i)*data.params.l(j);

data.freq = 0.5;

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


t = [0:0.1:10]';

for variable = 1:length(t)
    
    fprintf('%d\n',variable);
    
   trajd = Flat2state(testTraj(t(variable)),data.params); 
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

xd = [trajd.xL(n).x;trajd.xL(n).dx{1}];
xd = [xd;reshape(trajd.R, 9,1);trajd.Omega];
for i = 1:n
    xd = [xd;trajd.q(i).q];
end
for i = 1:n
    xd = [xd;vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
end
    
desiredstate(variable,:) = xd';
    
end


