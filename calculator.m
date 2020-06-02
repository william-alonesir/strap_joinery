%This script was designed for plotting the M-theta curve of the strap
%joinery. 
% Author: Zhuoyang Xin

% Material Property 
  E_per=1000;                                           %Mpa
  E_pl=8000;                                            %Mpa                           
% Sample Geometry
  I=2.355*10^7;                                         %mm^4
  L=600;                                                %mm
  h=200;                                                %mm
  a1=38.5;                                              %mm
  a2=161.5;                                             %mm
  w=43;                                                 %mm
  t=10;                                                 %mm

% 1st stage: linear elastic stage
  [Pc,angle_initial,k1,strain_initial]=linear_force(h,a1,a2,w,t,E_per,L);
  
% 2nd stage. Considering the geometry nonlinearity
  % Different cases for contact area (2flange1web,1flange1web,1flange only)
  [angle_c11,angle_c21,angle_c12,angle_c22]=critical_angle(angle_initial,h,a1,a2,w,t,E_per,L);
  
  % Relationship between moment and rotation
  [moment_load,moment_reaction,moment_strap,angle_matrix] =...,
      moment_rotation(angle_initial,angle_c11,angle_c21,angle_c12,angle_c22,h,a1,a2,w,t,E_per,L);
  
  % Moment-rotation Curve
  moment_load=[[0,k1*angle_initial],moment_load];
  angle_matrix=[[0,angle_initial],angle_matrix];
  figure(1)
  hold on
  plot(angle_matrix,moment_load)
  title('Moment-Rotation Theoretical Curve');
  xlabel ('Rotation Angle / rad')
  ylabel ('Moment / Nm')
  xlim([0,0.4]);
  ylim([0,1600]);
  hold off;
