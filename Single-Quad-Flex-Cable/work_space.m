% 
% syms th1 th2 l m g real;
% 
% fx = m*l*th1^2;
% 
% gx = th2/g;
% 
% 
% m_list_x = {'th1', 'x(1)' ; ...
%             'th2', 'x(2)' } ;
% m_list_params = {'l', 'params(1)' ; ...
%                  'm', 'params(2)' ; ...
%                  'g', 'params(3)'} ;
% 
% write_fcn_m('fcn_f_g.m',{'x', 'params'},[m_list_x;m_list_params],{fx,'fx'; gx,'gx'});
% 

% load('params.mat');
% traj = get_flat_traj(0);
% 
% [trajd] = get_desired_traj(traj,params);
% 
% time = [0:0.1:5]';
% 
% load3 = [];
% load2 = [];
% load1 = [];
% loadQ = [];
% 
% for i = 1:length(time)
%     [trajd(i)] = get_desired_traj(get_flat_traj(time(i)),params);
%     load3 = [load3;trajd(i).xL(3).x'];
%     load2 = [load2;trajd(i).xL(2).x'];
%     load1 = [load1;trajd(i).xL(1).x'];
%     loadQ = [loadQ;trajd(i).xQ.x'];
%     xQd(i,:) = trajd(i).xQ.x';
%     xLd(i,:) = 
%     qd(i,:) = 
%     Rd
% end
% 
% figure(1);
% plot3(load3(:,1),load3(:,2),load3(:,3));hold on;
% plot3(load2(:,1),load2(:,2),load2(:,3));hold on;
% plot3(load1(:,1),load1(:,2),load1(:,3));hold on;
% plot3(loadQ(:,1),loadQ(:,2),loadQ(:,3));hold on;
% axis equal;
% grid on;
% hold off;

% syms t real;
% t = 0;
% for t = 0:0.1:1
%     B = 2*[cos(t);0;sin(t)];
%     dB{1} = 2*[-sin(t);0;cos(t)];
%     dB{2} = 2*[-cos(t);0;-sin(t)];
%     dB{3} = 2*[sin(t);0;-cos(t)];
%     dB{4} = 2*[cos(t);0;sin(t)];
%     dB{5} = 2*[-sin(t);0;cos(t)];
%     dB{6} = 2*[-cos(t);0;-sin(t)];
%     dB{7} = 2*[sin(t);0;-cos(t)];
%     dB{8} = 2*[cos(t);0;sin(t)];
%     dB{9} = 2*[-sin(t);0;cos(t)];
%     dB{10} = 2*[-cos(t);0;-sin(t)];
% 
%     q2 = -[cos(t);0;sin(t)];
%     dq2{1} = -[-sin(t);0;cos(t)];
%     dq2{2} = -[-cos(t);0;-sin(t)];
%     dq2{3} = -[sin(t);0;-cos(t)];
%     dq2{4} = -[cos(t);0;sin(t)];
%     dq2{5} = -[-sin(t);0;cos(t)];
%     dq2{6} = -[-cos(t);0;-sin(t)];
%     dq2{7} = -[sin(t);0;-cos(t)];
%     dq2{8} = -[cos(t);0;sin(t)];
%     dq2{9} = -[-sin(t);0;cos(t)];
%     dq2{10} = -[-cos(t);0;-sin(t)];
% 
%     q= -B/norm(B);
%     dq = get_q_derivaties(B,norm(B),dB{:});
% 
%     compare = {q,dq{:};q2,dq2{:}};
%     
% end


%%
syms t B nB B(t) nB(t) real
dq = -B(t)/nB(t);

for i = 1:14
    dq = subs(diff(dq,t),diff(nB(t),t),dot(B(t),diff(B(t),t))/nB(t))
end
















