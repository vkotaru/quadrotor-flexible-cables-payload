% Minimum Snap Trajectory Generation 
% http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5980409
% 
% Author: vkotaru@andrew.cmu.edu
% Date: 7-June-2016
% Last Updated: 7-June-2016
% ================================================================
% Initialize Workspace:
clear;
clc;

% ======================
% Parameters & Constants
mu_r = 1;
mu_psi = 1;

n = 4; % degree of the polynomial
m = 1; % no. of segments
p = 4; % Number of elements

h1 = zeros(n+1);
h2(end) = 24^2;


h2 = zeros(n+1);
h2(end) = 24^2;

h3 = zeros(n+1);
h3(end) = 24^2;

h4 = zeros(n+1);
h4(3:end,3:end) = [4 6 8;6 12 18;8 18 144/5];

H = blkdiag(h1,h2,h3,h4);
c = zeros(4*m*(n+1),1);
f = zeros(size(c));

Aeq = [1 , zeros(1,19); zeros(1,5), 1, zeros(1,14);zeros(1,10),1,zeros(1,9);zeros(1,15),1,zeros(1,4);
    0,1,zeros(1,18);zeros(1,5),0,1,zeros(1,13); zeros(1,10),0,1,zeros(1,8);zeros(1,15),0,1,zeros(1,3);
    ones(1,5),zeros(1,15);zeros(1,5),ones(1,5),zeros(1,10); zeros(1,10),ones(1,5),zeros(1,5);zeros(1,15),ones(1,5);
    0,1,2,3,4,zeros(1,15);zeros(1,5),0:4,zeros(1,10);zeros(1,10),0:4,zeros(1,5);zeros(1,15),0:4];
Beq = [zeros(4,1);zeros(4,1);2;2;2;10;zeros(4,1)]; 

opts = optimoptions('quadprog',...
    'Algorithm','interior-point-convex','Display','iter');

[x fval eflag output lambda] = quadprog(H,f,[],[],Aeq,Beq,[],[],[],opts);



