function[traj] = get_flat_traj(t,varargin)
% Function to generate Desired Differentially-Flat trajectory 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 18-May-2016
% =====================================================================



traj.x=[2.*cos(0.15708E1.*t);0.25E1.*sin(0.15708E1.*t);0];
traj.dx{1}=[(-0.314159E1).*sin(0.15708E1.*t);0.392699E1.*cos(0.15708E1.*t);0];


traj.dx{2}=[(-0.49348E1).*cos(0.15708E1.*t);(-0.61685E1).*sin(0.15708E1.*t);0];


traj.dx{3}=[0.775157E1.*sin(0.15708E1.*t);(-0.968946E1).*cos(0.15708E1.*t);0];


traj.dx{4}=[0.121761E2.*cos(0.15708E1.*t);0.152202E2.*sin(0.15708E1.*t);0];


traj.dx{5}=[(-0.191262E2).*sin(0.15708E1.*t);0.239078E2.*cos(0.15708E1.*t);0];


traj.dx{6}=[(-0.300434E2).*cos(0.15708E1.*t);(-0.375543E2).*sin(0.15708E1.*t);0];


traj.dx{7}=[0.471921E2.*sin(0.15708E1.*t);(-0.589901E2).*cos(0.15708E1.*t);0];


traj.dx{8}=[0.741291E2.*cos(0.15708E1.*t);0.926614E2.*sin(0.15708E1.*t);0];


traj.dx{9}=[(-0.116442E3).*sin(0.15708E1.*t);0.145552E3.*cos(0.15708E1.*t);0];


traj.dx{10}=[(-0.182906E3).*cos(0.15708E1.*t);(-0.228633E3).*sin(0.15708E1.*t);0];


% traj.x = flip([0.1*t^2;0.2*t;0]);
% traj.dx{1} = flip([0.2*t;0.2;0]);
% traj.dx{2} = flip([0.2;0;0]);
% % traj.x = zeros(3,1);
% for i = 3:10
%    traj.dx{i} = zeros(3,1); 
% end

% traj.x = [0.001*t^5;0.001*t^4;0.001*t^3];
% traj.x = [0.005*t^4;0.004*t^3;0.003*t^2];
% traj.dx{1} = [0.02*t^3;0.012*t^2;0.006*t];
% traj.dx{2} = [0.06*t^2;0.024*t;0.006];
% traj.dx{3} = [0.12*t;0.024;0];
% traj.dx{4} = [0.12;0;0];
% for i = 5:10
%    traj.dx{i} = zeros(3,1); 
% end

% traj.x = flip([0.02*t^3;0.012*t^2;0.006*t]);
% traj.dx{1} = flip([0.06*t^2;0.024*t;0.006]);
% traj.dx{2} = flip([0.12*t;0.024;0]);
% traj.dx{3} = flip([0.12;0;0]);
% for i = 4:10
%    traj.dx{i} = zeros(3,1); 
% end

% traj.x = [0.02*t^3;0.012*t^2;0.006*t];
% traj.dx{1} = [0.06*t^2;0.024*t;0.006];
% traj.dx{2} = [0.12*t;0.024;0];
% traj.dx{3} = [0.12;0;0];
% for i = 4:10
%    traj.dx{i} = zeros(3,1); 
% end

end
