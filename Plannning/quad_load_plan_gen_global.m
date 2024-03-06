%% quad_load_plan_traj_global.m
% by Guofan Wu (gwu)
% -----------------------------------------------------------------------
% generate a trajectory offline for the single quadrotor with suspended load to 
% track from the start position, passing through a set of internal points and 
% return back to the start position x0 with the same yaw angle 
% Compute the coefficients based on a global QP 
% -----------------------------------------------------------------------
% This version still use the minimum snap trajectory generation 
% -----------------------------------------------------------------------
% quad_load_plan_gen_global(x0, ti, xi, delta_t,m, psi) 
%
% Input Argument: 
% x0                    initial start position with zero initial velocity 
% psi                   the constant yaw angle for these quadrotors
% delta_t               the sample time for trajectory generation 
% ti                    a list of internal time 
% xi                    a list of internal position that has to be reached at
% the specific internal point
% m                     the minimum degree of the polynomial for the load position, 
% should be at least 6 
% 
% 
% Return: store a matrix called "quad_load_plan.mat" for Simulink block to
% utilize in position control of the quadrotor
% All units here are in SI 
% Need to be debugged carefully for the first version 
% -----------------------------------------------------------------------

function quad_load_plan_gen_global(x0, ti, xi, delta_t, m, psi)
clc; 
close all; 

%% Preprocesssing 
% the input argument should be 
if nargin < 4
    error('Not enough input argument!'); 
end

if nargin == 4
    % default values for the polynomial degree and initial yaw angle 
    m = 6;
    psi = 0;  
end

if nargin == 5
    psi = 0; 
end

% number of the internal states  
n = length(ti); 
if size(xi, 2) ~= n
    error('The input argument is not consistent!');
end

if m < 6 
    disp('The input degree is not large enough! Use 6 instead now.');
    m = 6; 
end


%% Compute the coefficients using a QP 
disp('Generate the trajectory of the load position now.'); 

tol = 1e-6; % tolerance for debugging purpose 

% the time duration for each period  
tend = ti(end) + max(diff([0, ti]));
tt = [0, ti, tend]; 
duration = diff(tt); 
xx = [x0, xi]; 

% the degree of continuity that needs to be maintained 
k = 6; 
% the number of coefficients
p = k+1;    

% set up the boundary conditions for the overall QP 
optoption = optimset('Display', 'off', 'TolX', 1e-9);
qp_dim = (n+1) * (m+1);           % number of coefficients
cons_dim = n * (p+1) + 6;          % number of equality constraints: continuity at internal points 
H = eye(qp_dim); 
Aeq = zeros(cons_dim, qp_dim); 
Ip = eye(p); 
bxeq = zeros(cons_dim, 1); 
byeq = zeros(cons_dim, 1); 
bzeq = zeros(cons_dim, 1); 

%% add in the position constraints 
bxeq(1:n+1) = xx(1,:)';
byeq(1:n+1) = xx(2,:)';
bzeq(1:n+1) = xx(3,:)'; 

i = 1; 
for j = 1 : n+1
    Aeq(i, (m+1) * (j-1) + 1) = 1;
    i = i+1; 
end

%% add in the continuity constraints at internal states
for j = 1 : n+1
    [H1, Aeq1] = comp_H(m, duration(j), k); 
    
    % assembly the coefficients 
    index1 = (m+1) * (j - 1) + 5;
    index2 = (m+1) * j; 
    index3 = (m+1) * (j - 1) + 1; 
    index4 = (m+1) * j + p; 
    H(index1:index2, index1:index2) = H1; 
    
    if j < n+1
        Aeq(i:i+k, index3:index2) = Aeq1; 
        Aeq(i:i+k, index2+1:index4) = -Ip;
        i = i + p; 
    else
        % update the boundary on the end point
        bxeq(i:i+2) = [x0(1); 0; 0];
        byeq(i:i+2) = [x0(2); 0; 0];
        bzeq(i:i+2) = [x0(3); 0; 0]; 
        Aeq(i:i+2, index3:index2) = Aeq1(1:3, :);
        i = i + 3; 
    end
end

%% add in the boundary condition on start point: zero velocity and acceleration 
Aeq(i, 2) = 1;
i = i+1;
Aeq(i, 3) = 1;
clear i; 

ax = quadprog(H, [], [], [], Aeq, bxeq, [], [], [], optoption); 
ay = quadprog(H, [], [], [], Aeq, byeq, [], [], [], optoption); 
az = quadprog(H, [], [], [], Aeq, bzeq, [], [], [], optoption); 

ax = reshape(ax, m+1, n+1); 
ay = reshape(ay, m+1, n+1); 
az = reshape(az, m+1, n+1); 

%% interpolate between each iteration to get the corresponding values 
t_even = [];
pos_even = []; 
for j = 1 : n+1
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
    
    % compute the 2nd derivative of acceleration 
    d2acc_coe = zeros(length(tj), m+1); 
    for i = 4 : m
        d2acc_coe(:, i+1) = i * t_power(:, i-3) * (i-1) * (i-2) * (i-3); 
    end
    d2accj = [(d2acc_coe * ax(:, j))';
              (d2acc_coe * ay(:, j))'; 
              (d2acc_coe * az(:, j))'];
          
    % compute the 3rd derivative 
    d3acc_coe = zeros(length(tj), m+1); 
    for i = 5 : m
        d3acc_coe(:, i+1) = i * t_power(:, i-3) * (i-1) * (i-2) * (i-3) * (i-4); 
    end
    d3accj = [(d3acc_coe * ax(:, j))';
                    (d3acc_coe * ay(:, j))'; 
                    (d3acc_coe * az(:, j))'];
    
    % compute the 4th derivative 
    d4acc_coe = zeros(length(tj), m+1); 
    for i = 6 : m
        d4acc_coe(:, i+1) = i * t_power(:, i-3) * (i-1) * (i-2) * (i-3) * (i-4) * (i-5); 
    end
    d4accj = [(d4acc_coe * ax(:, j))';
                    (d4acc_coe * ay(:, j))'; 
                    (d4acc_coe * az(:, j))'];
    
    % concatenate the current position and velocity 
    if j == 1
        pos_even = posj;
        vel_even = velj;
        acc_even = accj;
        dacc_even = daccj;
        d2acc_even = d2accj; 
        d3acc_even = d3accj;
        d4acc_even = d4accj; 
        t_even = tj; 
    else
        % concatenate the rest of durations
        pos_even =   [pos_even, posj(:, 2:end)];
        vel_even =   [vel_even, velj(:, 2:end)];
        acc_even =   [acc_even, accj(:, 2:end)];
        dacc_even =  [dacc_even, daccj(:, 2:end)];
        d2acc_even = [d2acc_even, d2accj(:, 2:end)];
        d3acc_even = [d3acc_even, d3accj(:, 2:end)];
        d4acc_even = [d4acc_even, d4accj(:, 2:end)];
        t_even = [t_even; tj(2:end) + ti(j-1)]; 
    end
end

%% convert to the corresponding state of quadrotor
disp('The coefficients have been successfully generated.');
disp('Converting to the corresponding state using differential flatness.'); 

% mass properties of this quadrotor // update this if you are using another
% quadrotor
sys_prms.mQ = 0.545;
sys_prms.J = diag([0.00360481639,  0.00372100553,  0.00703101500]); 
sys_prms.g = 9.83; 
sys_prms.e3 = [0;0;1];
sys_prms.e1 = [1;0;0]; 

% mass properties of the load   // update this if you are using a new load
sys_prms.mL = 0.098;        % in kilogram 
sys_prms.L = 0.6;               % in meter

% allocate memory 
time_num = length(t_even); 
quad_ref = zeros(time_num, 15); 
load_ref = zeros(time_num, 9); 
R = zeros(3, 3, time_num); 
for j = 1 : time_num
    flat_output.psi = psi;
    flat_output.dpsi = 0;
    flat_output.d2psi = 0; 
    flat_output.x = pos_even(:, j); 
    flat_output.v = vel_even(:, j); 
    flat_output.acc = acc_even(:, j); 
    flat_output.dacc = dacc_even(:, j); 
    flat_output.d2acc = d2acc_even(:, j); 
    flat_output.d3acc = d3acc_even(:, j); 
    flat_output.d4acc = d4acc_even(:, j); 
    traj_val = Flat2Space(flat_output, sys_prms); 
    R(:, :, j) = traj_val.R; 
    quad_ref(j, :) = [traj_val.xQ', traj_val.vQ', traj_val.aQ', traj_val.euler_angle]; 
    load_ref(j, :) = [pos_even(:, j)', vel_even(:, j)', acc_even(:, j)']; 
end

quad_traj = [t_even, quad_ref]; 
load_traj = [t_even, load_ref]; 
save('quad_load_traj.mat', 'quad_traj', 'load_traj'); 

prompt = 'Do you want to test out the trajectory? (1: yes) > ';
test_flag = input(prompt); 

if test_flag == 1

    % debug this planner by plotting out the trajectory 
    fig_size = [1, 1, 6, 5]; 
    disp('Plotting out the trajectory now...');
    figure('Color', [1 1 1], 'Units', 'inches', 'Position', fig_size,...
           'PaperPositionMode', 'manual',   'PaperPosition', fig_size);
    
    pos_max = max(pos_even, [], 2); 
    pos_min = min(pos_even, [], 2); 
    diff_pos = pos_max - pos_min;
    
    axis_max = 1.2 * diff_pos + pos_min + 0.25;
    axis_min = -1.2 * diff_pos + pos_max - 0.25; 
    
    frame_rate = floor(1/delta_t);
    speed = 2.0; 
    voutput = VideoWriter('quad_plan_animation.avi');
    set(voutput, 'FrameRate',floor(frame_rate/speed));
    open(voutput); 
    
    
    quad_lw = 2.0;
    quad_clr = rand(1, 3); 
    for j = 1 : time_num      
        plot3(quad_traj(1:j, 2), quad_traj(1:j, 3), quad_traj(1:j, 4), 'r', 'LineWidth', 1.5); 
        hold on;
        plot3(quad_traj(j, 2), quad_traj(j, 3), quad_traj(j, 4), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', [1 0 0]); 
        draw_single_quad(quad_traj(j, 2:4)', R(:, :, j), quad_clr, quad_lw); 
        plot3(load_traj(j, 2), load_traj(j, 3), load_traj(j, 4), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0 0 1]); 
        plot3(load_traj(1:j, 2), load_traj(1:j, 3), load_traj(1:j, 4), 'b--', 'LineWidth', 2.5); 
        hline = line([load_traj(j, 2), quad_traj(j, 2)], [load_traj(j, 3), quad_traj(j, 3)], [load_traj(j, 4), quad_traj(j, 4)], 'Color', [0 0 0]);
        set(hline, 'LineWidth', 3.5);
        
        hold off;
        xlabel('X axis', 'FontSize', 16);
        zlabel('Y axis', 'FontSize', 16);
        ylabel('Z axis', 'FontSize', 16);
        box on; 
        set(gca, 'XLim', [axis_min(1), axis_max(1)], 'YLim', [axis_min(2), axis_max(2)], 'ZLim', [axis_min(3), axis_max(3)]); 
        drawnow;
        pause(delta_t); 
        cframe = getframe; 
        writeVideo(voutput, cframe); 
    end
    close(voutput); 
end

end

%% Compute the coefficients H for planning optimization 
% m is the number of polynomial degree, dt is the time duration and 
% k is the number of continuity that needs to preserved at intermediate
% point 
function [H, Aeq] = comp_H(m, dt, k) 
H = zeros(m-3); 
Aeq = zeros(k+1, m+1); 

% minimum snap 
for i = 4 : m
    for j = i : m
        alpha1 = factorial(i)/factorial(i-4);
        alpha2 = factorial(j)/factorial(j-4); 
        b = i + j - 7; 
        H(i-3, j-3) = alpha1 * alpha2 * dt^b/b; 
        H(j-3, i-3) = H(i-3, j-3); 
    end
end

% the degree of continuity that needs to be preserved 
for i = 0 : k
    alpha2 = factorial(i);
    for j = i : m
        b = j-i;
        alpha1 = factorial(j)/factorial(b); 
        Aeq(i+1, j+1) = alpha1 * dt^(b)/alpha2; 
    end
end


end

%% Inline function for converting the flat output to the corresond state
% for this case, the flat output is the corresponding load trajectory 
% and the yaw angle of load 
function [traj_val] = Flat2Space(flat_output, sys_params)
    mQ = sys_params.mQ ;
    mL = sys_params.mL ;
    J = sys_params.J ;
    L = sys_params.L ;
    g = sys_params.g ;
    e1 = sys_params.e1 ;
    e3 = sys_params.e3 ;
    
    xLd = flat_output.x;
    vLd = flat_output.v;
    aLd = flat_output.acc;
    daLd = flat_output.dacc ;
    d2aLd = flat_output.d2acc;
    d3aLd = flat_output.d3acc;
    d4aLd = flat_output.d4acc;
    
    Tp = -mL*(aLd + g*e3) ;
    norm_Tp = norm(Tp) ;
    p = Tp / norm_Tp ;
    
    dTp = -mL*daLd ;
    dnorm_Tp = 1/norm_Tp * (Tp' * dTp) ;
%     dp = (dTp*norm_Tp - Tp*dnorm_Tp)/norm_Tp^2 ;
%     dp = (dTp*norm_Tp - p*norm_Tp*dnorm_Tp)/norm_Tp^2 ;
    dp = (dTp - p*dnorm_Tp) / norm_Tp ;

    
    d2Tp = -mL*d2aLd ;
%     d2norm_Tp = ( (dot(dTp, dTp)+dot(Tp, d2Tp))*norm_Tp - dot(Tp, dTp)*dnorm_Tp ) / norm_Tp^2 ;
%     d2norm_Tp = ( (dot(dTp, dTp)+dot(Tp, d2Tp))*norm_Tp - dnorm_Tp*norm_Tp*dnorm_Tp ) / norm_Tp^2 ;
    d2norm_Tp = ( (dTp' * dTp) + (Tp' * d2Tp) - dnorm_Tp^2 ) / norm_Tp ;
%     d2p = ( (d2Tp*norm_Tp+dTp*dnorm_Tp - dTp*dnorm_Tp-Tp*d2norm_Tp)*norm_Tp^2 - (dp*norm_Tp^2)*2*norm_Tp*dnorm_Tp ) / norm_Tp^4 ;
%     d2p = ( (d2Tp - dp*dnorm_Tp-p*d2norm_Tp)*norm_Tp - (dp*norm_Tp)*2*norm_Tp*dnorm_Tp ) / norm_Tp^2 ;
    d2p = ( d2Tp - dp*dnorm_Tp - p*d2norm_Tp - dp*dnorm_Tp) / norm_Tp ;
    
    d3Tp = -mL*d3aLd ;
%     d3norm_Tp = ( (2*dot(dTp, d2Tp) + dot(dTp, d2Tp)+dot(Tp, d3Tp) - 2*dnorm_Tp*d2norm_Tp )*norm_Tp - d2norm_Tp*norm_Tp*dnorm_Tp  ) ...
%         / norm_Tp^2 ;
    d3norm_Tp = ( 2*dot(d2Tp, dTp) + dot(dTp, d2Tp)+dot(Tp, d3Tp) - 3*dnorm_Tp*d2norm_Tp) / norm_Tp ;
%     d3p = ( (d3Tp - d2p*dnorm_Tp-dp*d2norm_Tp - dp*d2norm_Tp-p*d3norm_Tp - 2*d2p*norm_Tp*dnorm_Tp-2*dp*dnorm_Tp^2-2*dp*norm_Tp*d2norm_Tp)*norm_Tp ...
%     - d2p*norm_Tp*dnorm_Tp) / norm_Tp^2 ;
    d3p = (d3Tp - d2p*dnorm_Tp-dp*d2norm_Tp - dp*d2norm_Tp-p*d3norm_Tp - d2p*dnorm_Tp-dp*d2norm_Tp - d2p*dnorm_Tp) / norm_Tp ;
    
    d4Tp = -mL*d4aLd ;
    d4norm_Tp = ( 2*dot(d3Tp, dTp)+2*dot(d2Tp, d2Tp) + dot(d2Tp, d2Tp)+dot(dTp, d3Tp) + dot(dTp, d3Tp)+dot(Tp, d4Tp) - 3*d2norm_Tp^2-3*dnorm_Tp*d3norm_Tp ...
        - d3norm_Tp*dnorm_Tp) / norm_Tp ;
    d4p = ( d4Tp - d3p*dnorm_Tp-d2p*d2norm_Tp - d2p*d2norm_Tp-dp*d3norm_Tp - d2p*d2norm_Tp-dp*d3norm_Tp - dp*d3norm_Tp-p*d4norm_Tp ...
        - d3p*dnorm_Tp-d2p*d2norm_Tp - d2p*d2norm_Tp-dp*d3norm_Tp - d3p*dnorm_Tp-d2p*d2norm_Tp - d3p*dnorm_Tp ) / norm_Tp ;
    
    axQ = aLd - L*d2p ;
    daxQ = daLd - L*d3p ;
    d2axQ = d2aLd - L*d4p ;
    
    b1d = e1 ;
    db1d = zeros(3,1) ;
    d2b1d = zeros(3,1) ;
    
    fb3 = mQ*(axQ+g*e3) - Tp ;
    norm_fb3 = norm(fb3) ;
    f = norm_fb3 ;
    b3 = fb3 / norm_fb3 ;
    b3_b1d = cross(b3, b1d) ;
    norm_b3_b1d = norm(b3_b1d) ;
    b1 = - cross(b3, b3_b1d) / norm_b3_b1d ;
    b2 = cross(b3, b1) ;
    R = [b1 b2 b3] ;
    
    dfb3 = mQ*(daxQ) - dTp ;
    dnorm_fb3 = dot(fb3, dfb3) / norm_fb3 ;
    db3 = (dfb3*norm_fb3 - fb3*dnorm_fb3) / norm_fb3^2 ;
    db3_b1d = cross(db3, b1d) + cross(b3, db1d) ;
    dnorm_b3_b1d = dot(b3_b1d, db3_b1d) / norm_b3_b1d ;
    db1 = (-cross(db3,b3_b1d)-cross(b3,db3_b1d) - b1*dnorm_b3_b1d) / norm_b3_b1d ;
    db2 = cross(db3, b1) + cross(b3, db1) ;
    dR = [db1 db2 db3] ;
    Omega = vee(R'*dR) ;
    
    d2fb3 = mQ*(d2axQ) - d2Tp ;
    d2norm_fb3 = (dot(dfb3, dfb3)+dot(fb3, d2fb3) - dnorm_fb3*dnorm_fb3) / norm_fb3 ;
    d2b3 = ( (d2fb3*norm_fb3+dfb3*dnorm_fb3 - dfb3*dnorm_fb3-fb3*d2norm_fb3)*norm_fb3^2 - db3*norm_fb3^2*2*norm_fb3*dnorm_fb3 ) / norm_fb3^4 ;
    d2b3_b1d = cross(d2b3, b1d)+cross(db3, db1d) + cross(db3, db1d)+cross(b3, d2b1d) ;
    d2norm_b3_b1d = ( (dot(db3_b1d,db3_b1d)+dot(b3_b1d,d2b3_b1d))*norm_b3_b1d - dot(b3_b1d, db3_b1d)*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
    d2b1 = ( (-cross(d2b3,b3_b1d)-cross(db3,db3_b1d) - cross(db3,db3_b1d)-cross(b3,d2b3_b1d) - db1*dnorm_b3_b1d-b1*d2norm_b3_b1d )*norm_b3_b1d - db1*norm_b3_b1d*dnorm_b3_b1d ) / norm_b3_b1d^2 ;
    d2b2 = cross(d2b3, b1)+cross(db3, db1) + cross(db3, db1)+cross(b3, d2b1) ;
    d2R = [d2b1 d2b2 d2b3] ;
    dOmega = vee( dR'*dR + R'*d2R ) ; %vee( dR'*dR + R'*d2R, true ) ;
    
    M = J*dOmega + cross(Omega, J*Omega) ;
    
    % compute the velocity and acceleration of the corresponding quadrotor 
    xQ = xLd - L*p ;
    vQ = vLd - L * dp; 
    aQ = aLd - L * d2p; 
    
    % compute the corresponding Euler angles 
    [phi, theta, psi] = rotm2euler(R); 

    % compute the corresponding angular velocity 
    ang_vel = [1, sin(phi) * tan(theta), cos(phi) * tan(theta);
               0,            cos(phi),   -sin(phi);
               0,            sin(phi)*sec(theta), cos(phi)*sec(theta)] * Omega; 
    
    traj_val.xQ = xQ;
    traj_val.vQ = vQ;
    traj_val.aQ = aQ;
    traj_val.q = p;
    traj_val.dq =dp;
    traj_val.d2q =d2p;
    traj_val.R = R;
    traj_val.Omega =Omega;
    traj_val.dOmega = dOmega;
    traj_val.f = f;
    traj_val.M = M;
    traj_val.euler_angle = [ang_vel; M]'; 
    traj_val.phi = phi;
    traj_val.theta = theta;
    traj_val.psi = psi; 
end

%% Convert rotation matrix to Euler angles 
function [phi, theta, psi] = rotm2euler(R)
if abs(R(3, 1)) ~= 1
    theta = -asin(R(3, 1)); 
    psi = atan2(R(2,1)/cos(theta), R(1,1)/cos(theta));
    phi = atan2(R(3,2)/cos(theta), R(3,3)/cos(theta));
else if R(3, 1) == -1
        theta = -pi/2;
        psi = 0;
        phi = atan2(R(1, 2), R(1,3));
    else
        theta = pi/2; 
        psi = 0;
        phi = atan2(-R(1,2), -R(1,3)); 
    end
end
end


