function[traj] =testTraj(t,plannedTraj)
% Function to generate Desired Differentially-Flat trajectory 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 9-Sep-2016
% =====================================================================

% Parameters
path = plannedTraj.path;
time = plannedTraj.time;



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
