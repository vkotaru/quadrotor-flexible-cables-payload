%% Function to  study the nature of 'Quadrotor with flexible Cable system'
% with varying frequencies. 
% 
%
% Author: vkotaru@andrew.cmu.edu
% Date: 2-June-2016

% Last Updated: 14-Jun-2016

%% INITIALZING WORKSPACE
% =====================
% profile on
% Clear workspace
% clear; 
close all; 
% clc;

% Add Paths
% Geometric Control Toolbox
addpath('../GeoControl-Toolbox/');

%% System constants and parameters
% ===============================
data.params.mQ = 0.5 ;
data.params.J = diag([0.557, 0.557, 1.05]*10e-2);
data.params.g = 9.81 ;
data.params.e1 = [1;0;0] ;
data.params.e2 = [0;1;0] ;
data.params.e3 = [0;0;1] ;

data.params.m = [0.05,0.05,0.05,0.05,0.5]; % 0.1*ones(1,5); % mass of each link
data.params.l = 0.25*ones(1,5); % length of each link 
data.params.n = 5; % No. of links suspended

% data.params.m = [0.14]; % 0.1*ones(1,5); % mass of each link
% data.params.l = 1.25*ones(1,1); % length of each link 
% data.params.n = 1; % No. of links suspended


n = data.params.n;
data.params.M00 = data.params.mQ + sum(data.params.m(1:data.params.n));
data.params.M0i = @(i) sum(data.params.m(i:data.params.n))*data.params.l(i); 
data.params.Mi0 = data.params.M0i;
data.params.Mij = @(i,j) sum(data.params.m(max(i,j):data.params.n))*data.params.l(i)*data.params.l(j);

frequency = 0.1:0.1:1;

%%

export = cell(length(frequency),2);

for freq = 1:length(frequency)
    
    for j = 1:2
        
        if j == 1
            data.params.m = [0.05,0.05,0.05,0.05,0.05]; % 0.1*ones(1,5); % mass of each link
            data.params.l = 0.25*ones(1,5); % length of each link 
            data.params.n = 5; % No. of links suspended
        elseif j == 2
%             keyboard;
            data.params.m = [0.25]; % 0.1*ones(1,5); % mass of each link
            data.params.l = 1.25*ones(1,1); % length of each link 
            data.params.n = 1; % No. of links suspended
        end
    
        xL = [];
        xQ = [];
        R = [];
        f = [];
        M = [];
        time = [];
        l = [];

        for t = 0:0.1:1

        fprintf('Frequency %.2f,n = %d, time %.2f\n',frequency(freq),data.params.n,t);

            [trajd] = get_desired_traj(get_flat_circ_traj(t,frequency(freq)),data.params);

            xL_ = trajd.xL(data.params.n).x ;
            vL_ = trajd.xL(data.params.n).dx{1};
            R_ = trajd.R;
            Omega_ = trajd.Omega;

            q = [];
            dq = [];
    %         omega = [];
            for i = 1:data.params.n
                q = [q,trajd.q(i).q];
                dq = [dq, trajd.q(i).dq{1}];
    %             omega = [omega, vec_cross(trajd.q(i).q,trajd.q(i).dq{1})];
            end

            xQ_ = xL_ - sum(repmat(data.params.l,3,1).*q,2);
            vQ_ = vL_ - sum(repmat(data.params.l,3,1).*dq,2);
            l_ = norm(xQ_-xL_);
            
            xL = [xL; xL_',vL_'];
            R = [R; reshape(R_,1,9),Omega_'];
            xQ = [xQ; xQ_',vQ_'];
            f = [f;trajd.f];
            M = [M;trajd.M'];
            l = [l;l_];
%             time = [time;t];

        end

        output.xL = xL;
        output.xQ = xQ;
        output.R = R;
        output.f = f;
        output.M = M;
        output.t = t;
        output.data = data;
        output.l = l;
        
        export{freq,j} = output;
        
    end
end

%%
save('export_freq3.mat','export');
run('temp');

%% IGNORE 
% % close all;
% % freq = [1;0.5;0.25;0.125];
% % 
% % 
% % for i = 1:length(freq)
% %    
% %     output{i} =  main_quad_flexible(freq(i));
% %     
% % end
% % 
% % save('outputn1.mat','output');
% % 
% % 
% % figure;
% % for i = 1:length(output)
% % %     subplot(2,2,i);
% %     plot(output{i}.t,output{i}.f); hold on;
% %     grid on; 
% %     
% %     
% % end
% %   legend('1','.5','.25','.125');
% 
% clear;
% close all;
% files = dir('output/*.mat');
% list = {'rs','go','b+','y*','yx'};
% list2 = {};
% 
% figure;
% for i = 1:length(files)
%     
%     clear('output');
%     load(fullfile('output',files(i).name));
% %     
% %     f(1,i) = max(output{1}.f);
% %     f(2,i) = max(output{2}.f);
% %     f(3,i) = max(output{3}.f);
% %     f(4,i) = max(output{4}.f);
% %     
% %     freq(1,i) = output{1}.freq;
% %     freq(2,i) = output{2}.freq;
% %     freq(3,i) = output{3}.freq;
% %     freq(4,i) = output{4}.freq;
%     
% %     plot(t,f); hold on; grid on;
% 
%     clear('x');clear('xQ');clear('xL');
%     x = output{1}.x;
%     l = output{1}.data.params.l';
%     n = output{1}.data.params.n;    
%     for j = 1:size(x,1)
% 
%         xL(j,:) = x(j,1:3);
%         q = reshape(x(j,19:3*n+18),3,n)';
% 
%         xQ(j,:) = xL(j,:) - sum(repmat(l(1:n),1,3).*q,1);
% 
%     end
% %     plot3(xL(:,1),xL(:,2),xL(:,3),'-r','LineWidth',1);hold on; grid on;
%     plot3(xQ(:,1),xQ(:,2),xQ(:,3),'-g');hold on; grid on;
%     
%     clear('x');clear('xQ');clear('xL');
%     x = output{2}.x;
%     l = output{2}.data.params.l';
%     n = output{2}.data.params.n;    
%     for j = 1:size(x,1)
% 
%         xL(j,:) = x(j,1:3);
%         q = reshape(x(j,19:3*n+18),3,n)';
% 
%         xQ(j,:) = xL(j,:) - sum(repmat(l(1:n),1,3).*q,1);
% 
%     end
% %     plot3(xL(:,1),xL(:,2),xL(:,3),'-r','LineWidth',1);hold on; grid on;
%     plot3(xQ(:,1),xQ(:,2),xQ(:,3),'-b');hold on; grid on;
%     
%     clear('x');clear('xQ');clear('xL');
%     x = output{3}.x;
%     l = output{3}.data.params.l';
%     n = output{3}.data.params.n;    
%     for j = 1:size(x,1)
% 
%         xL(j,:) = x(j,1:3);
%         q = reshape(x(j,19:3*n+18),3,n)';
% 
%         xQ(j,:) = xL(j,:) - sum(repmat(l(1:n),1,3).*q,1);
% 
%     end
% %     plot3(xL(:,1),xL(:,2),xL(:,3),'-r','LineWidth',1);hold on; grid on;
%     plot3(xQ(:,1),xQ(:,2),xQ(:,3),'-y');hold on; grid on;
%     
%     clear('x');clear('xQ');clear('xL');
%     x = output{4}.x;
%     l = output{4}.data.params.l';
%     n = output{4}.data.params.n;    
%     for j = 1:size(x,1)
% 
%         xL(j,:) = x(j,1:3);
%         q = reshape(x(j,19:3*n+18),3,n)';
% 
%         xQ(j,:) = xL(j,:) - sum(repmat(l(1:n),1,3).*q,1);
% 
%     end
%     plot3(xL(:,1),xL(:,2),xL(:,3),'-r','LineWidth',1);hold on; grid on;
%     plot3(xQ(:,1),xQ(:,2),xQ(:,3),'-m');hold on; grid on; axis equal;
%  
% %     n(i) = output{1}.data.params.n;
% %     list2{i} = num2str(n(i));
% %     
% %     scatter(freq(:,i),f(:,i),list{i});grid on;hold on;   
% end
% 
% % legend(list2);

