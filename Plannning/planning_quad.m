function[varargout] = planning_quad(ti,xi,n,r,initialGuessScale)
%% planning
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

%% quadrotor parameters
% Quadrotor parameters
data.params.mQ = 0.5 ;
data.params.J = diag([0.557, 0.557, 1.05]*10e-2);
data.params.g = 9.81 ;
data.params.e1 = [1;0;0] ;
data.params.e2 = [0;1;0] ;
data.params.e3 = [0;0;1] ;

% data.params.m = [0.5]; % mass of each link
% data.params.l = 0.5*ones(1,1); % length of each link 
% data.params.n = 1; % No. of links suspended

data.params.m = [0.25, 0.25]; % mass of each link
data.params.l = 0.25*ones(1,2); % length of each link 
data.params.n = 2; % 

data.params.M00 = data.params.mQ + sum(data.params.m(1:data.params.n));
data.params.M0i = @(i) sum(data.params.m(i:data.params.n))*data.params.l(i); 
data.params.Mi0 = data.params.M0i;
data.params.Mij = @(i,j) sum(data.params.m(max(i,j):data.params.n))*data.params.l(i)*data.params.l(j);
data.params.imax = 4 + 2*data.params.n;

if ~nargin
%    ti = [0,0.5,1];
%    xi = [0,0,0;2,0,2;0,2,0];
%    n = 6;
%    r = 4;
   ti = [0,1];
   xi = [0,0,0;2,0,0];
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
params.t_start = params.ti(1:end-1);
params.t_stop = params.ti(2:end);


data.params.xi = params.xi;
%% Objective Function and Constraints
Aeq = [];
Beq = [];
dBeq = [];

for i = 1:l
    
    h{i} = [];
    dA{i} = zeros(m+1,m*n1,r);
    dB = zeros((m+1)*r,1);
        
    for j = 1:m
        
        derv = poly_diff(N,r);
        coeff  = poly_int(derv,ti(j),ti(j+1),'coeff');
        coeff = [zeros(1,r),coeff]; %error here
        h{i} = blkdiag(h{i},coeff'*coeff);
        
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
Beq = [Beq;dBeq];

H = blkdiag(h{:});
x0 = initialGuessScale*ones(length(H),1);
params.H = H;

%% Optimization
disp('optimizing...');
% opts = optimoptions('quadprog',...
%     'Algorithm','interior-point-convex','Display','iter');
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

fun = @(x) obj_func(x,params,data); % Objective function
nonlcon =@(x) quad_constraints(x,params,data); % Non-linear constraints

% [x fval eflag output lambda] = quadprog(H,f,[],[],Aeq,Beq,[],[],[],opts);
[x,fval,exitflag,output]  = fmincon(fun,x0',[],[],Aeq,Beq,[],[],nonlcon,options);

%% Post processing
X = reshape(x,n1,m,l);
params.X = X;
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

%% OUTPUT
if nargout == 1
    varargout{1}.X = X;
    varargout{1}.params = params;
    varargout{1}.quad_data = data;
end

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

function[fvalue] = obj_func(x,params,data)

fvalue = params.mu_r*x*params.H*x';
% delT = [params.ti(1):0.01:params.ti(end)]';
% X = reshape(x,params.n1,params.m,params.l);
% 
% sum_f = 0;
% sum_M = 0;
% for j = 1:length(delT)
%     t = delT(j); 
%     [states] = map_t2state(t,X,params,data);
%     
%     sum_f = states.f^2 + sum_f;
%     sum_M = norm(states.M)^2 + sum_M;
% end
% 
% fvalue =  1*sum_f + 1*sum_M;

end

%% Non-linear constraints
function[c,ceq] = quad_constraints(x,params,data)
c = [];
% decision variables
X = reshape(x,params.n1,params.m,params.l);

% flat outputs n thier higher derivatives
t = params.ti(end);
time = (t*params.N).^params.powers;
path.x = [time*X(:,end,1);time*X(:,end,2);time*X(:,end,3)];

for i = 1:data.params.imax
    temp = (t*params.N(1:end-i)).^params.powers(1:end-i).*poly_diff(params.N,i);
    time = [zeros(1,i),temp];
    path.dx{i} = [time*X(:,end,1);time*X(:,end,2);time*X(:,end,3)];
end

% flatoutputs to state space
states = Flat2state(path,data.params);
% R = states.R - eye(3);
% final velocity and acceleartion of quadrotor are kept zero
ceq = [states.xQ.dx{1};states.xQ.dx{2}];%reshape(R,numel(R),1);states.Omega];


% delT = [params.ti(1):0.1:params.ti(end)]';
% s = size(delT,1);
% temp = ((delT*params.N(1:end-2)).^repmat(params.powers(1:end-2),s,1)).*repmat(poly_diff(params.N,2),s,1);
% time = [zeros(s,2),temp];
% 
% delT = [params.ti(1):0.1:params.ti(end)]';
% for j = 1:length(delT)
%     t = delT(j); 
% 
%     [states] = map_t2state(t,X,params,data);
%     c = [c;states.f-50;norm(states.M)-25];
% end
% 


end

function [state] = map_t2state(t,X,params,data)

iter = find((t>=params.t_start)&(t<params.t_stop));
if t == params.ti(end)
    iter = params.m;
end

time = (t*params.N).^params.powers;
path.x = [time*X(:,iter,1);time*X(:,iter,2);time*X(:,iter,3)];

for i = 1:data.params.imax
    temp = (t*params.N(1:end-i)).^params.powers(1:end-i).*poly_diff(params.N,i);
    time = [zeros(1,i),temp];
    path.dx{i} = [time*X(:,iter,1);time*X(:,iter,2);time*X(:,iter,3)];
end
state = Flat2state(path,data.params);

end




