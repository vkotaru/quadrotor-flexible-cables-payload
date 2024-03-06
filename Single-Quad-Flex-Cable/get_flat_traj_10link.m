function[traj] = get_flat_traj_10link(t,freq)
% Function to generate Desired Differentially-Flat trajectory 
% 
% Author: vkotaru@andrew.cmu.edu
% Last Updated: 18-May-2016
% =====================================================================



traj.x=[0.2E1.*cos(freq.*pi.*t);0.25E1.*sin(freq.*pi.*t);0];
traj.dx{1}=[(-0.628319E1).*freq.*sin(freq.*pi.*t);0.785398E1.*freq.*cos(freq.*pi.* ...
  t);0];


traj.dx{2}=[(-0.197392E2).*freq.^2.*cos(freq.*pi.*t);(-0.24674E2).*freq.^2.*sin( ...
  freq.*pi.*t);0];


traj.dx{3}=[0.620126E2.*freq.^3.*sin(freq.*pi.*t);(-0.775157E2).*freq.^3.*cos( ...
  freq.*pi.*t);0];


traj.dx{4}=[0.194818E3.*freq.^4.*cos(freq.*pi.*t);0.243523E3.*freq.^4.*sin(freq.* ...
  pi.*t);0];


traj.dx{5}=[(-0.612039E3).*freq.^5.*sin(freq.*pi.*t);0.765049E3.*freq.^5.*cos( ...
  freq.*pi.*t);0];


traj.dx{6}=[(-0.192278E4).*freq.^6.*cos(freq.*pi.*t);(-0.240347E4).*freq.^6.*sin( ...
  freq.*pi.*t);0];


traj.dx{7}=[0.604059E4.*freq.^7.*sin(freq.*pi.*t);(-0.755073E4).*freq.^7.*cos( ...
  freq.*pi.*t);0];


traj.dx{8}=[0.189771E5.*freq.^8.*cos(freq.*pi.*t);0.237213E5.*freq.^8.*sin(freq.* ...
  pi.*t);0];


traj.dx{9}=[(-0.596182E5).*freq.^9.*sin(freq.*pi.*t);0.745227E5.*freq.^9.*cos( ...
  freq.*pi.*t);0];


traj.dx{10}=[(-0.187296E6).*freq.^10.*cos(freq.*pi.*t);(-0.23412E6).*freq.^10.*sin( ...
  freq.*pi.*t);0];


traj.dx{11}=[0.588408E6.*freq.^11.*sin(freq.*pi.*t);(-0.73551E6).*freq.^11.*cos( ...
  freq.*pi.*t);0];


traj.dx{12}=[0.184854E7.*freq.^12.*cos(freq.*pi.*t);0.231067E7.*freq.^12.*sin(freq.* ...
  pi.*t);0];


traj.dx{13}=[(-0.580735E7).*freq.^13.*sin(freq.*pi.*t);0.725919E7.*freq.^13.*cos( ...
  freq.*pi.*t);0];


traj.dx{14}=[(-0.182443E8).*freq.^14.*cos(freq.*pi.*t);(-0.228054E8).*freq.^14.*sin( ...
  freq.*pi.*t);0];




end
