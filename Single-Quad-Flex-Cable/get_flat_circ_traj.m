function[traj] =get_flat_circ_traj(t,freq)
% Function to generate Desired Differentially-Flat trajectory 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 18-May-2016
% =====================================================================



traj.x=[0.5E0.*cos(2.*freq.*pi.*t);0.5E0.*sin(2.*freq.*pi.*t);0];
traj.dx{1}=[(-0.314159E1).*freq.*sin(2.*freq.*pi.*t);0.314159E1.*freq.*cos(2.* ...
  freq.*pi.*t);0];


traj.dx{2}=[(-0.197392E2).*freq.^2.*cos(2.*freq.*pi.*t);(-0.197392E2).*freq.^2.* ...
  sin(2.*freq.*pi.*t);0];


traj.dx{3}=[0.124025E3.*freq.^3.*sin(2.*freq.*pi.*t);(-0.124025E3).*freq.^3.*cos( ...
  2.*freq.*pi.*t);0];


traj.dx{4}=[0.779273E3.*freq.^4.*cos(2.*freq.*pi.*t);0.779273E3.*freq.^4.*sin(2.* ...
  freq.*pi.*t);0];


traj.dx{5}=[(-0.489631E4).*freq.^5.*sin(2.*freq.*pi.*t);0.489631E4.*freq.^5.*cos( ...
  2.*freq.*pi.*t);0];


traj.dx{6}=[(-0.307645E5).*freq.^6.*cos(2.*freq.*pi.*t);(-0.307645E5).*freq.^6.* ...
  sin(2.*freq.*pi.*t);0];


traj.dx{7}=[0.193299E6.*freq.^7.*sin(2.*freq.*pi.*t);(-0.193299E6).*freq.^7.*cos( ...
  2.*freq.*pi.*t);0];


traj.dx{8}=[0.121453E7.*freq.^8.*cos(2.*freq.*pi.*t);0.121453E7.*freq.^8.*sin(2.* ...
  freq.*pi.*t);0];


traj.dx{9}=[(-0.763113E7).*freq.^9.*sin(2.*freq.*pi.*t);0.763113E7.*freq.^9.*cos( ...
  2.*freq.*pi.*t);0];


traj.dx{10}=[(-0.479478E8).*freq.^10.*cos(2.*freq.*pi.*t);(-0.479478E8).*freq.^10.* ...
  sin(2.*freq.*pi.*t);0];


traj.dx{11}=[0.301265E9.*freq.^11.*sin(2.*freq.*pi.*t);(-0.301265E9).*freq.^11.*cos( ...
  2.*freq.*pi.*t);0];


traj.dx{12}=[0.18929E10.*freq.^12.*cos(2.*freq.*pi.*t);0.18929E10.*freq.^12.*sin(2.* ...
  freq.*pi.*t);0];


traj.dx{13}=[(-0.118935E11).*freq.^13.*sin(2.*freq.*pi.*t);0.118935E11.*freq.^13.* ...
  cos(2.*freq.*pi.*t);0];


traj.dx{14}=[(-0.747288E11).*freq.^14.*cos(2.*freq.*pi.*t);(-0.747288E11).* ...
  freq.^14.*sin(2.*freq.*pi.*t);0];




end
