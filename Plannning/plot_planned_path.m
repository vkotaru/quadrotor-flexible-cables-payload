function plot_planned_path(input)
%
%
%% Quadrotor parameters
data = input.quad_data;
% flat outputs n thier higher derivatives
imax = 2*data.params.n + 4;


%% path planning
params = input.params;

X = input.X;
powers = params.powers;
N = params.N;
ti = params.ti;
m = params.m;

t_start = ti(1:end-1);
t_stop = ti(2:end);

total_time = ti(1):0.1:ti(end);
points = [];


for j = 1:length(total_time) 
    fprintf('iteration %d of %d\n',j,length(total_time));
    t = total_time(j);
    time = (t*params.N).^params.powers;

    iter = find((t>=t_start)&(t<t_stop));
    if t == ti(end)
        disp('THE END');
        iter = m;
    end
    
    
    path.x = [time*X(:,iter,1);time*X(:,iter,2);time*X(:,iter,3)]; 

    for i = 1:imax
        temp = (t*params.N(1:end-i)).^params.powers(1:end-i).*poly_diff(params.N,i);
        time = [zeros(1,i),temp];
        path.dx{i} = [time*X(:,iter,1);time*X(:,iter,2);time*X(:,iter,3)]; 
    end
    
    trajd0 = Flat2state(path,data.params);
    xL0 = trajd0.xL(data.params.n).x ;
    vL0 = trajd0.xL(data.params.n).dx{1};
    aL0 = trajd0.xL(data.params.n).dx{2};
    R0 = trajd0.R;
    Omega0 = trajd0.Omega;

    q0 = [];
    dq0 = [];
    d2q0 = [];
    omega0 = [];
    for i = 1:data.params.n
        q0 = [q0,trajd0.q(i).q];
        dq0 = [dq0, trajd0.q(i).dq{1}];
        d2q0 = [d2q0, trajd0.q(i).dq{2}];
        omega0 = [omega0, vec_cross(trajd0.q(i).q,trajd0.q(i).dq{1})];
    end
    xQ0 = xL0 - sum(repmat(data.params.l,3,1).*q0,2);
    vQ0 = vL0 - sum(repmat(data.params.l,3,1).*dq0,2);
    aQ0 = vL0 - sum(repmat(data.params.l,3,1).*d2q0,2);
    
    xQd(j,:) = [xQ0',vQ0',aQ0'];
    xLd(j,:) = [xL0',vL0',aL0'];
    xd(j,:) = [xL0; vL0; reshape(R0,9,1); Omega0; reshape(q0,numel(q0),1); reshape(omega0,numel(omega0),1)];

    f(j,1) = trajd0.f;
    M(j,:) = trajd0.M';
end
% 

%% Plots
keyboard;
t = total_time';
figure;
subplot(2,2,1);
plot(t,xQd(:,1),'r',t,xQd(:,2),'g',t,xQd(:,3),'b');grid on;
legend('x','y','z');
title('Quad position');
subplot(2,2,2);
plot(t,xQd(:,4),'r',t,xQd(:,5),'g',t,xQd(:,6),'b');grid on;
legend('x','y','z');
title('Quad Velocity');
subplot(2,2,3);
plot(t,xQd(:,7),'r',t,xQd(:,8),'g',t,xQd(:,9),'b');grid on;
legend('x','y','z');
title('Quad Acceleration');

figure;
subplot(2,2,1);
plot(t,xLd(:,1),'r',t,xLd(:,2),'g',t,xLd(:,3),'b');grid on;
title('Load position');
subplot(2,2,2);
plot(t,xLd(:,4),'r',t,xLd(:,5),'g',t,xLd(:,6),'b');grid on;
title('Load Velocity');
subplot(2,2,3);
plot(t,xLd(:,7),'r',t,xLd(:,8),'g',t,xLd(:,9),'b');grid on;
title('Load Acceleration');

figure;
subplot(1,2,1);
plot(t,f,'LineWidth',2);grid on;
title('thrust');
subplot(1,2,2);
plot(t,M(:,1),'r',t,M(:,2),'g',t,M(:,3),'b','LineWidth',1);grid on;
title('moment');


animate_3dquad_flexcable(t, xd, data.params,0);











end