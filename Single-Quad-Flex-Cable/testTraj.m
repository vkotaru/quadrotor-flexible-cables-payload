function[traj] =testTraj(t)
% Function to generate Desired Differentially-Flat trajectory 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 18-May-2016
% =====================================================================

% Parameters
f1 = 1/4; f2 = 1/5; f3 = 1/7;
ax = 2; ay = 2.5; az = 1.5;



traj.x=[t;t^2;t^3];
traj.dx{1}=[1; 2*t; 3*t^2];


traj.dx{2}=[0; 2; 6*t];


traj.dx{3}=[0;0;6];


traj.dx{4}=[0;0;0];


traj.dx{5}=[0;0;0];


traj.dx{6}=[0;0;0];


traj.dx{7}=[0;0;0];


traj.dx{8}=[0;0;0];


traj.dx{9}=[0;0;0];


traj.dx{10}=[0;0;0];


traj.dx{11}=[0;0;0]; 

traj.dx{12}=[0;0;0];

traj.dx{13}=[0;0;0];


traj.dx{14}=[0;0;0];




end
