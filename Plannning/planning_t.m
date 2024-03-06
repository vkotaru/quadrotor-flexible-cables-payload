function planning_t(ti,xi,n,r)
%% planning THIS CODE IS WRONG 
% 
%
% ---------------------------------------------------------------------
% planning for quadrotor with load suspended via 'n' link flexible cable.
% Initial and and final conditions specified.
% ---------------------------------------------------------------------
% Input Arguments:
%   ti [1 x (m+1)] - (m>1)time intevals from 't0' (initial) to 'tf' (final time)
%   xi [(m+1) x 3] - points at different time instances
%   n - degree of the polynomials
%   r - order of the derivative to be optimized. For Minimum snap r=4
% ---------------------------------------------------------------------
% Author: Prasath Kotaru (vkotaru@andrew.cmu.edu)
% Date: June-25-2016
% Last Updated: July-13-2016
% ======================================================================

%% INITIALIZING WORKSPACE
% =======================
% Clear workspace
% clear; 
close all; 
% clc;

% Add Paths
% Geometric Control Toolbox
addpath('../GeoControl-Toolbox/');

if ~nargin
   ti = [0,0.5,1];
   xi = [0,0,0;2,0,2;0,2,0];
   n = 6;
   r = 4;
end

%% INITIALIZING PARAMETERS
% ========================
mu_r = .1;
l = size(xi,2);     % [x,y,z] or [x,y,z,psi]
m = length(ti)-1;   % Total no. of polynominals btn the internals i.e., to--t1--..--tf  -- is a polynomials
n1 = n+1;           % No. of coefficients in each polynomials i.e., degree of polynomial + constant
%% 
h = cell(1,l);
N = ones(1,n+1);
powers = 0:n;

%% Parameters
% ---------
params.mu_r = mu_r;
params.l = l;
params.n = n;
params.n1 = n1;
params.N = N;
params.powers = powers;
params.ti = ti;
params.xi = xi;
params.r = r;
params.m = m;

%%
Aeq = [];
Beq = [];
dBeq = [];

for i = 1:l
    
    dA{i} = zeros(m+1,m*n1,r);
    dB = zeros((m+1)*r,1);    
    for j = 1:m
        % position constraints
        Aeq = blkdiag(Aeq,[(ti(j)*N).^powers;(ti(j+1)*N).^powers]);
        Beq = [Beq;xi(j,i);xi(j+1,i)];
        
        % derivative constraints
        for k = 1:r % Upto r(th) derivative equality 
            derv_1 = poly_diff(N,k);
            powers_1 = powers(1:end-k);
            N_1 = N(1:end-k);
            
            dA{i}(j,n1*(j-1)+1:n1*j,k) = -[zeros(1,k) (ti(j)*N_1).^powers_1.*derv_1];
            dA{i}(j+1,n1*(j-1)+1:n1*j,k) = [zeros(1,k) (ti(j+1)*N_1).^powers_1.*derv_1];
        end       
    end
    
    % restructuring the derivative constraints
    dA_{i} = [];
    for kk = 1:r
       dA_{i} = [dA_{i};dA{i}(:,:,kk)];  
    end
    dBeq = [dBeq;dB];   
end

dAeq = blkdiag(dA_{:});

Aeq = [Aeq;dAeq];
Aeq = [Aeq, zeros(size(Aeq,1),1)];
Beq = [Beq;dBeq];

% H = blkdiag(h{:});
f = ones(size(Aeq,2),1);

% opts = optimoptions('quadprog',...
%     'Algorithm','interior-point-convex','Display','iter');
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

fun = @(x) obj_func(x,params);

% [x fval eflag output lambda] = quadprog(H,f,[],[],Aeq,Beq,[],[],[],opts);
[x,fval,exitflag,output]  = fmincon(fun,f',[],[],Aeq,Beq,[],[],[],options);

X = reshape(x(1:end-1),n1,m,l);
params.X = X;
ti = ti*x(end); params.ti = ti;
t = ti(1):0.01:ti(end);
points = [];
for i = 1:length(t) 
    points = [points;construct_path(t(i),params)];
end
% 
figure;
plot3(points(:,1),points(:,2),points(:,3),'b','LineWidth',1);hold on;
scatter3(xi(:,1),xi(:,2),xi(:,3),'sr');
grid on;


end

function [path_point] = construct_path(t,params)
X = params.X;
powers = params.powers;
N = params.N;
ti = params.ti;
m = params.m;

if t < min(ti) || t > max(ti)
    error('t is not in the range of ti');
end

time = (t*N).^powers;

for i = 1:m
    if t >= ti(i) && t<ti(i+1)
        path_point = [time*X(:,i,1),time*X(:,i,2),time*X(:,i,3)]; 
    elseif t == ti(end)
        disp('THE END');
        path_point = [time*X(:,m,1),time*X(:,m,2),time*X(:,m,3)]; 
    end
end


end

%%
function[fvalue] = obj_func(x,params)
% parameters
powers = params.powers;
N = params.N;
ti = params.ti;
m = params.m;
l = params.l;
r = params.r;
mu_r = params.mu_r;

alpha = x(end);
x = x(1:end-1);
ti = ti/alpha;

for i = 1:l
    
    h{i} = [];
    for j = 1:m
        
        derv = poly_diff(N,r);
        coeff  = poly_int(derv,ti(j),ti(j+1),'coeff');
        coeff = [zeros(1,r),coeff]; %error here
        h{i} = blkdiag(h{i},coeff'*coeff);
           
    end

end

H = blkdiag(h{:});

fvalue = abs(alpha^(2*r-1))*mu_r*x*H*x';

end


