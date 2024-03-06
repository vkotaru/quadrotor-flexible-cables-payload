function[traj] = get_flat_traj(t,varargin)
% Function to generate Desired Differentially-Flat trajectory 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 18-May-2016
% =====================================================================



traj.x=[2.*(1+(-1).*cos((1/2).*pi.*t));0.25E1.*sin((2/5).*pi.*t);0.15E1.*cos(( ...
  2/7).*pi.*t)];
traj.dx{1}=[pi.*sin((1/2).*pi.*t);0.314159E1.*cos((2/5).*pi.*t);(-0.13464E1).*sin(( ...
  2/7).*pi.*t)];


traj.dx{2}=[(1/2).*pi.^2.*cos((1/2).*pi.*t);(-0.394784E1).*sin((2/5).*pi.*t);( ...
  -0.120852E1).*cos((2/7).*pi.*t)];


traj.dx{3}=[(-1/4).*pi.^3.*sin((1/2).*pi.*t);(-0.4961E1).*cos((2/5).*pi.*t); ...
  0.108477E1.*sin((2/7).*pi.*t)];


traj.dx{4}=[(-1/8).*pi.^4.*cos((1/2).*pi.*t);0.623418E1.*sin((2/5).*pi.*t); ...
  0.973685E0.*cos((2/7).*pi.*t)];


traj.dx{5}=[(1/16).*pi.^5.*sin((1/2).*pi.*t);0.78341E1.*cos((2/5).*pi.*t);( ...
  -0.873978E0).*sin((2/7).*pi.*t)];


traj.dx{6}=[(1/32).*pi.^6.*cos((1/2).*pi.*t);(-0.984463E1).*sin((2/5).*pi.*t);( ...
  -0.784481E0).*cos((2/7).*pi.*t)];


traj.dx{7}=[(-1/64).*pi.^7.*sin((1/2).*pi.*t);(-0.123711E2).*cos((2/5).*pi.*t); ...
  0.704148E0.*sin((2/7).*pi.*t)];


traj.dx{8}=[(-1/128).*pi.^8.*cos((1/2).*pi.*t);0.15546E2.*sin((2/5).*pi.*t); ...
  0.632042E0.*cos((2/7).*pi.*t)];


traj.dx{9}=[(1/256).*pi.^9.*sin((1/2).*pi.*t);0.195357E2.*cos((2/5).*pi.*t);( ...
  -0.56732E0).*sin((2/7).*pi.*t)];


traj.dx{10}=[(1/512).*pi.^10.*cos((1/2).*pi.*t);(-0.245493E2).*sin((2/5).*pi.*t);( ...
  -0.509225E0).*cos((2/7).*pi.*t)];




end
