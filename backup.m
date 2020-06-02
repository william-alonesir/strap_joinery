% Material Property 
  E_per=1000;
  E_pl=8000;%Elastic Modulus
  
% Sample Geometry
  I=2.355*10^7; %Inertia
  L=600;% Length 
  % Cross section Dimension
  h=200;
  a1=38.5;
  a2=161.5;
  w=43;
  t=10;
%   
% 1st stage. 0<P<Pe 
  % Obtaining critical value  
    % obtaining strain value
    syms strain
    fun_strain = (w*a1^2+ t*(a2-a1)*(a1+a2)+w*a1*(a2+h))*strain*E_per/(2*h)==12000*0.9;
    s=double(solve(fun_strain,strain));
    % obtaining Pe
    fun1 = @(x) w*x*s*E_per.*(h-x)/h;
    fun2 = @(x) t*x*s*E_per.*(h-x)/h;
    Pe = (6000*0.9*h-integral(fun1,0,a1)-integral(fun2,a1,a2)-integral(fun1,a2,h))/L;
    % Flexural deflection at loading point
    % v=-(Pe*L^3/(6*E_pl*I)-L^2*Pe*L/(2*E_pl*I));
    % Total Deflection at loading point
    d=s*h*L/h;
    % Rotation Angle Measured from Loading Point 
    % angle=rad2deg(L^2*Pe/(2*E*I)+asin(s*h*L/h/L));
    angle_rad=asin(d/L);
    angle_deg=rad2deg(angle_rad);
    % Stiffness in 1st stage
    k1=Pe*L/ angle_rad/1e3;
    disp(['The initial rotational stiffness is ',num2str(k1),' Nm/rad']);
  % Obtaining stiffness for both conditions
  % ke=Pe/s;
  % kf=Pe/v;
    
% 2nd Stage
  %The case between the start to seperate and zero contact area between
  %beam_1 and beam_2
  % Bottom Strap Elongation vs Rotation
  fs1=@(x)400*(sqrt(3+cos(x)+2*sin(x))-2);        
  x = linspace(0,0.6);
  figure (1);
  hold on;
  title('Bottom Strap Extension-Rotation','FontSize',12);
  xlabel('Rotation Angle /rad');
  ylabel('Extension /mm');
  plot(x, fs1(x))
  hold off;
  % Top Strap Elongation vs Rotation
  fs2=@(x)800-200*(sqrt(10+6*cos(x)));                     % in deg
  figure (2);
  hold on;
  title('Top Strap Extension-Rotation','FontSize',12);
  xlabel('Rotation Angle /rad');
  ylabel('Extension /mm');
  plot(x, fs2(x))
  hold off;
  % Strap Constitutive Curve
  fs=@(ex)-7*10^(-8)*ex.^4 + 2*10^(-5)*ex.^3 - 0.0017*ex.^2 + 0.0989*ex; %
  ex = linspace(0,140);
  figure (3);
  hold on;
  plot(ex, fs(ex))
  hold off;
  % Relationship between Rotation and maximum strain 
  alpha1=@(x)atan((L*sin(x)-h*cos(x)+h)/(L*cos(x)+h*sin(x)+h));     % in rad
  alpha2=@(x)atan((L*sin(x)/(L*cos(x)+h)));                         % in rad
  
  syms strain2 strain3
%   x2=atan(strain2*h/a2);
%   fun_strain2 = (t*(a2-a1)^2+w*a1*(2*a2-a1))*strain2*E_per/(2*a2)==2*fs(48.7+fs1(rad2deg(x2)))*cos(alpha1(x2))+2*fs(48.7-fs2(rad2deg(x2)))*cos(alpha2(x2));
    fun_strain2 = (t*(a2-a1)^2+w*a1*(2*a2-a1))*strain2*E_per/(2*a2)==12000*0.9;
    s2=double(solve(fun_strain2,strain2));
    x2=atan(s2*h/a2);
    fun_strain3 = w*a1*strain3*E_per/2==12000*0.9;
    s3=double(solve(fun_strain3,strain3));
    x3=atan(s3*h/a1);
%   integral(fun2,a1,a2)
  fs111=400*(sqrt(3+cos(deg2rad(3))+2*sin(deg2rad(3)))-2);
 
%   E_per*s2*A*b/2=2*fs(48.7+fs1(x))*cos(alpha1)+2*fs(48.7-fs2(x))*cos(alpha2)
  %x=atan(s2*h/b);                                       %The relationship between rotation and s2
  moment_matrix=[];
  angle_matrix=[];
  moment_reaction=[];
  for x=angle_rad:0.00002:x2
      syms len
      fun2_strain = (w*(len-a2)^2+t*(a2-a1)*(len-a2+len-a1)+w*a1*(2*len-a1))*(len*tan(x)/h)*E_per./(2*len)==12000*0.9;      
      len4=double(solve(fun2_strain,len));
      strain4=len4(2)*tan(x)/h;
      fun3 = @(locate) w*locate*strain4*E_per.*(len4(2)-locate)/len4(2);
      fun4 = @(locate) t*locate*strain4*E_per.*(len4(2)-locate)/len4(2);
      d1=sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(pi/4+x+atan(3))./(400*sqrt(3+cos(x)+2*sin(x)));
      d2=600*200*sin(pi-x)./(200*(sqrt(10+6*cos(x))));
      Mpc2=(2*fs(48.7+fs1(x)).*d1-2*fs(48.7-fs2(x)).*d2);
      M1=@(x)integral(fun3,0,a1)+integral(fun4,a1,a2)+integral(fun3,a2,len4(2));
      Mp2=Mpc2-M1(x)/1000;
      moment_matrix=[moment_matrix,Mp2];
      moment_reaction=[moment_reaction,M1(x)/1000];
      angle_matrix=[angle_matrix,x];
  end
  
% for x=x2:0.0002:0.03
for x=x2:0.0006:x3
      syms len2
      fun3_strain = (t*(len2-a1)^2+w*a1*(2*len2-a1))*(len2*tan(x)/h)*E_per./(2*len2)==12000*0.9;      
      len5=double(solve(fun3_strain,len2));
      strain5=len5(2)*tan(x)/h;
      fun5 = @(locate) w*locate*strain5*E_per.*(len5(2)-locate)/len5(2);
      fun6 = @(locate) t*locate*strain5*E_per.*(len5(2)-locate)/len5(2);
      d1=sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(pi/4+x+atan(3))./(400*sqrt(3+cos(x)+2*sin(x)));
      d2=600*200*sin(pi-x)./(200*(sqrt(10+6*cos(x))));
      Mpc2=(2*fs(48.7+fs1(x)).*d1-2*fs(48.7-fs2(x)).*d2);
      M1=@(x)integral(fun5,0,a1)+integral(fun6,a1,len5(2));
      Mp2=Mpc2-M1(x)/1000;
      moment_matrix=[moment_matrix,Mp2];
      moment_reaction=[moment_reaction,M1(x)/1000];
      angle_matrix=[angle_matrix,x];
end


for x=x3:0.005:0.4
    syms len3
    fun4_strain=w*len3*(len3*tan(x)/h)*E_per/2==12000*0.9;
    len6=double(solve(fun4_strain,len3));
      strain6=len6(2)*tan(x)/h;
      fun7 = @(locate) w*locate*strain6*E_per.*(len6(2)-locate)/len6(2);
      d1=sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(pi/4+x+atan(3))./(400*sqrt(3+cos(x)+2*sin(x)));
      d2=600*200*sin(pi-x)./(200*(sqrt(10+6*cos(x))));
      Mpc2=(2*fs(48.7+fs1(x)).*d1-2*fs(48.7-fs2(x)).*d2);
      M1=@(x)integral(fun7,0,len6(2));
      Mp2=Mpc2-M1(x)/1000;
      moment_matrix=[moment_matrix,Mp2];
      moment_reaction=[moment_reaction,M1(x)/1000];
      
      angle_matrix=[angle_matrix,x];
end
  % E_per*s2*A*b/2==2*fs(48.7+fs1(x))*cos(alpha1)+2*fs(48.7-fs2(x))*cos(alpha2);
  figure(4)b
  plot(angle_matrix,moment_matrix)
  figure (5)
  plot(angle_matrix,moment_reaction)
 
% % 3rd Stage
%   % Force vs Rotation angle
%   x=linspace(0,0.6);
%   d1=sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(pi/4+x+atan(3))./(400*sqrt(3+cos(x)+2*sin(x)));
%   d2=600*200*sin(pi-x)./(200*(sqrt(10+6*cos(x))));
%   dp=L*cos(x)+h*sin(x);
%   Mp4=(2*fs(48.7+fs1(x)).*d1-2*fs(48.7-fs2(x)).*d2);
%   % fp=@(t)600*(2*fs(48.7+400*(sqrt(3+cos(deg2rad(t))+2*sin(deg2rad(t)))-2)).*(sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(deg2rad(45+t+rad2deg(atan(3))))./(400*sqrt(3+cos(deg2rad(t))+2*sin(deg2rad(t)))))...,
%   %     -2*fs(48.7+200*(sqrt(10+6*cos(deg2rad(t))))-800).*(600*200*sin(deg2rad(180-t))./(200*(sqrt(10+6*cos(deg2rad(t)))))))...,
%   %     ./(L*cos(deg2rad(t))+h*sin(deg2rad(t)));
%   figure (5);
%   hold on;
%   plot (x,Mp4);
%   xlim([0 0.4])
%   ylim([0 1800])
%   title('Moment-Rotation Curve Stage 3','FontSize',12);
%   xlabel('Rotation Angle/rad');
%   ylabel ('Moment/Nm');
%   hold off;
  

% Conclusion
figure(6)
  hold on
  x=linspace(0,angle_rad);           % 1st stage
  plot(x,k1*x);                      
  plot(angle_matrix,moment_matrix);  % 2nd stage
  title('Moment-Rotation Theoretical Curve');
  xlabel('Rotation Angle / rad');
  ylabel('Moment / Nm');
  hold off


% % Find critical value
%   for x=angle_rad:0.00002:x2
%       syms len
%       fun2_strain = (w*(len-a2)^2+t*(a2-a1)*(len-a2+len-a1)+w*a1*(2*len-a1))*(len*tan(x)/h)*E_per./(2*len)==12000*0.9;      
%       len4=double(solve(fun2_strain,len));
%       strain4=len4(2)*tan(x)/h;
%       fun3 = @(locate) w*locate*strain4*E_per.*(len4(2)-locate)/len4(2);
%       fun4 = @(locate) t*locate*strain4*E_per.*(len4(2)-locate)/len4(2);
%       M1=@(x)integral(fun3,0,len4(2)-a2)+integral(fun4,len4(2)-a2,len4(2)-a1)+integral(fun3,len4(2)-a1,len4(2));
%       Mp2=(6000*0.9*h-M1(x))/1000;
%       d1=sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(deg2rad(45+x+rad2deg(atan(3))))./(400*sqrt(3+cos(deg2rad(x))+2*sin(deg2rad(x))));
%       d2=600*200*sin(deg2rad(180-x))./(200*(sqrt(10+6*cos(deg2rad(x)))));
%       dp=L*cos(deg2rad(x))+h*sin(deg2rad(x));
%       Mp4c=(2*fs(48.7+fs1(x)).*d1-2*fs(48.7-fs2(x)).*d2);
%       moment_matrixc=[moment_matrixc,Mp2-Mp4c];
%       angle_matrixc=[angle_matrixc,x];
%       if abs(Mp2-Mp4c)<0.55
%           xc=x;
%           break
%       end
%   end
  
%   figure(6)
%   hold on
%   x=linspace(0,angle_rad);           % 1st stage
%   plot(x,k1*x);                      
%   plot(angle_matrixc,moment_matrix(1:46));  % 2nd stage
%   x=linspace(rad2deg(xc),35);
%   plot(deg2rad(x),Mp4);
%   xlim([0,0.35]);
%   hold off
%   
%   
%   
% % Energy Check
% %   figure(6);
% %    hold on;
% %   angledisp=@(x)600*sin(deg2rad(x))+200-200*cos(deg2rad(x));
% %   plot(x,angledisp(x));
% %    plot(x,11*x);
% %    hold on;
% %   figure(5);
% %   hold on;
% %   angle=[];
% %   residualenergy=[];
% %   for theta=0:0.1:35
% %       %Up=@(x)fp(x).*(600*cos(deg2rad(x))+200*sin(deg2rad(x)));
% %       Up=@(x)fp(x)*11;
% %       dpp=600*sin(deg2rad(theta))+200-200*cos(deg2rad(theta));
% %       Uext=integral(Up,0, theta);
% %         % Us1=@(x)fs(400*(sqrt(3+cos(deg2rad(x))+2*sin(deg2rad(x)))-2)+48.7).*((2*cos(deg2rad(x))-sin(deg2rad(x)))*200./sqrt(3+cos(deg2rad(x))+2*sin(deg2rad(x)))+48.7);
% %       a2=48.7+400*(sqrt(3+cos(deg2rad( theta))+2*sin(deg2rad( theta)))-2);
% %       a22=800-200*(sqrt(10+6*cos(deg2rad( theta))));
% %       Us1output=2*integral(fs,48.7,a2);
% %       Us2output=2*integral(fs,48.7,48.7-a22);
% %       RU=Uext-Us1output-Us2output;
% %       plot(theta,RU);
% %       angle(length(angle)+1)=[theta];
% %       residualenergy(length(residualenergy)+1)=[RU];
% %       plot(angle,residualenergy);
% %   end
