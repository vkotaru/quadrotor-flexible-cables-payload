%% quad_plan_traj.m
% by Guofan Wu (gwu)
% -----------------------------------------------------------------------
% generate a trajectory offline for the single quadrotor to track from the
% start position, passing through a set of internal points and return back to the
% start position x0 through yaw angle 
% Propogate the coefficients from one interval to the next
% has the problem of overfitting 
% -----------------------------------------------------------------------
% quad_plan_traj(x0, delta_t, ti, xi, psi) 
%
% Input Argument: 
% x0                 initial start position with zero initial velocity 
% psi                the constant yaw angle for these quadrotors
% delta_t        the sample time for trajectory generation 
% ti                 a list of internal time 
% xi                 a list of internal position that has to be reached at
% the specific internal point
% m               x the maximum degree of the polynomial for the position, 
% should be at least 4  
% 
% Return: store a matrix called "quad_plan.mat" for Simulink block to
% utilize  
% The file "quad_plan.mat" is utilized in position control of the quadrotor
% Need to be debugged carefully for the first version 
% -----------------------------------------------------------------------

function quad_plan_gen(x0, delta_t, ti, xi, m, psi)
clc; 
close all; 

% tolerance for debugging
tol = 1e-6; 

% number of the internal states  
n = length(ti); 
if size(xi, 2) ~= n
    disp('The input argument is not consistent!');
    return; 
end
disp('Generate the trajectory now.'); 

% the time duration for each period  
tend = ti(end) + 4; % Why +4 here ?
tt = [0, ti, tend]; 
duration = diff(tt); 

% coefficients for the intermediate polynomials 
ax = zeros(m+1, n+1); 
ay = zeros(m+1, n+1); 
az = zeros(m+1, n+1); 

% set up the boundary condition at the initial condition  
ax(1, 1) = x0(1);
ay(1, 1) = x0(2);
az(1, 1) = x0(3); 

optoption = optimset('Display', 'off', 'TolX', 1e-8);

for j = 1 : n
    % get the current value of weight function H and constant vector f 
    [H, Aeq, Beq] = comp_H(m, duration(j)); 
    fx = comp_f(m, ax(5, j), duration(j)); 
    fy = comp_f(m, ay(5, j), duration(j)); 
    fz = comp_f(m, az(5, j), duration(j)); 
    
    Aeq1 = Aeq(6:m+1); 
    beqx = xi(1, j)-Aeq(1:5) * ax(1:5, j); 
    beqy = xi(2, j)-Aeq(1:5) * ay(1:5, j); 
    beqz = xi(3, j)-Aeq(1:5) * az(1:5, j); 
    
    % compute the corresponding coefficient 
    ux = quadprog(H, fx, [], [], Aeq1, beqx, [], [], [], optoption);  
    uy = quadprog(H, fy, [], [], Aeq1, beqy, [], [], [], optoption); 
    uz = quadprog(H, fz, [], [], Aeq1, beqz, [], [], [], optoption); 
    
    ax(6:end, j) = ux;
    ay(6:end, j) = uy;
    az(6:end, j) = uz;
    
    % propagate forward based on continuity 
    % position 
    ax(1, j+1) = xi(1, j);
    ay(1, j+1) = xi(2, j);
    az(1, j+1) = xi(3, j);
    
    % its higher order derivatives
    ax(2:5, j+1) = Beq * ax(:, j);
    ay(2:5, j+1) = Beq * ay(:, j);
    az(2:5, j+1) = Beq * az(:, j);
    
end

disp('The coefficients have been successfully generated.');
% Return to the original position with zero velocity 
j = n+1;
[H, Aeq, ~] = comp_H(m, duration(j)); 
fx = comp_f(m, ax(5, j), duration(j)); 
fy = comp_f(m, ay(5, j), duration(j)); 
fz = comp_f(m, az(5, j), duration(j)); 

Aeq1 = Aeq(6:m+1); 
beqx = x0(1)-Aeq(1:5) * ax(1:5, j); 
beqy = x0(2)-Aeq(1:5) * ay(1:5, j); 
beqz = x0(3)-Aeq(1:5) * az(1:5, j); 

ux = quadprog(H, fx, [], [], Aeq1, beqx, [], [], [], optoption);  
uy = quadprog(H, fy, [], [], Aeq1, beqy, [], [], [], optoption); 
uz = quadprog(H, fz, [], [], Aeq1, beqz, [], [], [], optoption); 

ax(6:end, j) = ux;
ay(6:end, j) = uy;
az(6:end, j) = uz;


% interpolate between each iteration to get the corresponding values 
t_even = [];
pos_even = []; 
for j = 1 : n
    tj = 0 : delta_t : duration(j); 
    tj = tj'; 
    
    % compute position and its higher order derivatives 
    t_power = zeros(length(tj), m+1); 
    for i = 1 : m+1
        t_power(:, i) = tj.^(i-1); 
    end
    posj = [(t_power * ax(:, j))';
            (t_power * ay(:, j))'; 
            (t_power * az(:, j))'];
    
    % compute the coefficients for velocity 
    vel_coe = zeros(length(tj), m+1); 
    for i = 1 : m
        vel_coe(:, i+1) = i * t_power(:, i); 
    end
    velj = [(vel_coe * ax(:, j))';
            (vel_coe * ay(:, j))'; 
            (vel_coe * az(:, j))'];
    
    % compute the acceleration
    acc_coe = zeros(length(tj), m+1); 
    for i = 2 : m
        acc_coe(:, i+1) = i * t_power(:, i-1) * (i-1); 
    end
    accj = [(acc_coe * ax(:, j))';
            (acc_coe * ay(:, j))'; 
            (acc_coe * az(:, j))'];
    
    % compute the derivative of acceleration
    dacc_coe = zeros(length(tj), m+1); 
    for i = 3 : m
        dacc_coe(:, i+1) = i * t_power(:, i-2) * (i-1) * (i-2); 
    end
    daccj = [(dacc_coe * ax(:, j))';
             (dacc_coe * ay(:, j))'; 
             (dacc_coe * az(:, j))'];
    
    % compute the 2nd derivative 
    d2acc_coe = zeros(length(tj), m+1); 
    for i = 4 : m
        acc_coe(:, i+1) = i * t_power(:, i-3) * (i-1) * (i-2) * (i-3); 
    end
    d2accj = [(d2acc_coe * ax(:, j))';
              (d2acc_coe * ay(:, j))'; 
              (d2acc_coe * az(:, j))'];
        
    % concatenate the current position and velocity 
    if j == 1
        pos_even = posj;
        vel_even = velj;
        acc_even = accj;
        dacc_even = daccj;
        d2acc_even = d2accj; 
        t_even = tj; 
    else
        % debug for correctness 
        if norm(pos_even(:, end) - posj(:, 1)) > tol
            disp('The distance is not continuous at the connecting point'); 
        end
        if norm(vel_even(:, end) - velj(:, 1)) > tol
            disp('The distance is not continuous at the connecting point'); 
        end
        if norm(acc_even(:, end) - accj(:, 1)) > tol
            disp('The distance is not continuous at the connecting point'); 
        end
        if norm(dacc_even(:, end) - daccj(:, 1)) > tol
            disp('The distance is not continuous at the connecting point'); 
        end
        if norm(d2acc_even(:, end) - d2accj(:, 1)) > tol
            disp('The distance is not continuous at the connecting point'); 
        end
        pos_even =   [pos_even, posj(:, 2:end)];
        vel_even =   [vel_even, velj(:, 2:end)];
        acc_even =   [acc_even, accj(:, 2:end)];
        dacc_even =  [dacc_even, daccj(:, 2:end)];
        d2acc_even = [d2acc_even, d2accj(:, 2:end)];
        t_even = [t_even; tj(2:end) + ti(j-1)]; 
    end
end


% debug this planner by plotting out the trajectory 
fig_size = [1, 1, 7, 6]; 
figure('Color', [1 1 1], 'Units', 'inches', 'Position', fig_size,...
       'PaperPositionMode', 'manual',   'PaperPosition', fig_size);
plot3(pos_even(1,:), pos_even(2,:), pos_even(3, :), 'r-', 'LineWidth', 2.3);
hold on;
plot3([xi(1,:), x0(1)], [xi(2,:), x0(2)], [xi(3,:), x0(3)], 'ro', 'MarkerFaceColor', [1 0 0], 'MarkerSize', 8); 
hold off; 
end

%% Compute the coefficients H for planning optimization 
function [H, Aeq, Beq] = comp_H(m, dt) 
H = zeros(m-4); 
Aeq = zeros(1, m+1); 
Beq = zeros(4, m+1); 

for i = 5 : m
    for j = i : m
        alpha1 = factorial(i)/factorial(i-4);
        alpha2 = factorial(j)/factorial(j-4); 
        b = i + j - 7; 
        H(i-4, j-4) = alpha1 * alpha2 * dt^b/b; 
        H(j-4, i-4) = H(i-4, j-4); 
    end
end

for j = 1 : m+1
    Aeq(j) = dt^(j-1); 
end

for i = 1 : 4
    alpha2 = factorial(i);
    for j = i : m
        b = j-i;
        alpha1 = factorial(j)/factorial(b); 
        Beq(i, j+1) = alpha1 * dt^(b)/alpha2; 
    end
end


end

%% Compute the constant term b using QP 
function [f] = comp_f(m, a4, dt) 
alpha1 = factorial(4); 
f = zeros(m-4, 1); 

for i = 5 : m
    alpha2 = factorial(i)/factorial(i-4);
    b = i - 3;
    f(i-4) = 2*alpha1 * alpha2 * a4 * dt^b/b; 
end
end
