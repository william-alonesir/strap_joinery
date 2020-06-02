% Plots experimental, numerical and theoretical plots in one figure

% moment_a1=xlsread('Testing Summary.xlsx','E3:E16356');
% % timber_a1=xlsread('Testing Summary.xlsx','F3:F16356');
% % strap_a1_ang=xlsread('Testing Summary.xlsx','C3:C16088');
% % strap_a1=xlsread('Testing Summary.xlsx','G3:G16088');
% angle_a1=xlsread('Testing Summary.xlsx','C3:C16356');
% moment_a2=xlsread('Testing Summary.xlsx','M3:M11014');
% % timber_a2=xlsread('Testing Summary.xlsx','N3:N11014');
% % strap_a2=xlsread('Testing Summary.xlsx','O3:O11014');
% angle_a2=xlsread('Testing Summary.xlsx','K3:K11014');
% moment_a3=xlsread('Testing Summary.xlsx','U3:U14061');
% % timber_a3=xlsread('Testing Summary.xlsx','V3:V14061');
% % strap_a3=xlsread('Testing Summary.xlsx','W3:W14061');
% angle_a3=xlsread('Testing Summary.xlsx','S3:S14061');
% % 
% moment_aba=xlsread('Testing Summary.xlsx','Abaqus_new','D3:D57');
% angle_aba=xlsread('Testing Summary.xlsx','Abaqus_new','C3:C57');

% Moment-rotation curve
figure (1)
hold on
plot(angle_matrix,moment_load,'DisplayName','Theoretical Results');
plot(angle_a1,moment_a1,'DisplayName','Experimental A1');
plot(angle_a2,moment_a2,'DisplayName','Experimental A2');
plot(angle_a3,moment_a3,'DisplayName','Experimental A3');
plot(angle_aba,moment_aba,'DisplayName','Numerical Results');
xlim([0,0.4]);
ylim([0,2100]);
% title('Numerical Curve of Moment-Rotation Relationship')
title('Strap Joinery Moment-Rotation Curve Comparison','fontsize',12)
xlabel('Rotaion Angle / rad');
ylabel('Moment / Nm');
legend
hold off

%Timber Rotation Curve
% figure (2)
% hold on
% plot(angle_a1,timber_a1);
% plot(angle_a2,timber_a2);
% plot(angle_a3,timber_a3);
% xlim([0,0.4]);
% ylim([-300,60]);
% title('Timber Strain-Rotation Relationship')
% xlabel('Rotaion Angle / rad');
% ylabel('Strain / micro');
% hold off

% Strap Rotation Curve
% figure (3)
% hold on
% plot(strap_a1_ang,strap_a1);
% plot(angle_a2,strap_a2);
% plot(angle_a3,strap_a3);
% fs1=@(x)400*(sqrt(3+cos(x)+2*sin(x))-2);  % plot theoretical strap_extension      
% x = linspace(0,0.4);
% plot(x, fs1(x)/800*1e6)
%  
% xlim([0,0.4]);
% ylim([0,80000]);
% title('Bottom Strain-Rotation Relationship')
% title('Comparison of Bottom Strap Strain-Rotation')
% xlabel('Rotaion Angle / rad');
% ylabel('Strain / micro');
% hold off

