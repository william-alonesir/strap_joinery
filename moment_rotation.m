function [moment_load,moment_reaction,moment_strap,angle_matrix] = moment_rotation(angle_initial,angle_c11,angle_c21,angle_c12,angle_c22,h,a1,a2,w,t,E_per,L)
% Obtaining relationship between moment and rotation in 2nd stage
%   The stage is seperated into three cases which different contact area
n=100;
alpha1=@(x)atan((L*sin(x)-h*cos(x)+h)/(L*cos(x)+h*sin(x)+h));       % in rad
alpha2=@(x)atan((L*sin(x)/(L*cos(x)+h)));                           % in rad
fs=@(x)1000*(2*fs_b(x)*cos(alpha1(x))+2*fs_t(x)*cos(alpha2(x)));    % in N
d_b=@(x)sqrt(L^2+h^2)*sqrt(h^2+h^2)*sin(pi/4+x+atan(3))./(400*sqrt(3+cos(x)+2*sin(x)));
d_t=@(x)600*200*sin(pi-x)./(200*(sqrt(10+6*cos(x))));
Mpc=@(x)(2*fs_b(x)*d_b(x)-2*fs_t(x)*d_t(x));
moment_load=[];
moment_reaction=[];
moment_strap=[];
angle_matrix=[];
% 1st case
xsi=angle_initial;
for k=1:n+1
      syms len
      fun1_strain = (w*(len-a2)^2+t*(a2-a1)*(len-a2+len-a1)+w*a1*(2*len-a1))*(len*tan(xsi)/h)*E_per./(2*len)==fs(xsi);      
      len1=double(solve(fun1_strain,len));
      strain1=len1(2)*tan(xsi)/h;
      fun1 = @(locate) w*locate*strain1*E_per.*(len1(2)-locate)/len1(2);
      fun2 = @(locate) t*locate*strain1*E_per.*(len1(2)-locate)/len1(2);
      M_r=@(x)(integral(fun1,0,a1)+integral(fun2,a1,a2)+integral(fun1,a2,len1(2)))/1000;
      Mp=Mpc(xsi)-M_r(xsi);
      moment_load=[moment_load,Mp];
      moment_reaction=[moment_reaction,M_r(xsi)];
      moment_strap=[moment_strap,Mpc(xsi)];
      angle_matrix=[angle_matrix,xsi];
      xsi=xsi+(angle_c11-angle_initial)/n;
end
% 2nd case
xsi=angle_c12;
for k=1:n+1
      syms len
      fun2_strain = (t*(len-a1)^2+w*a1*(2*len-a1))*(len*tan(xsi)/h)*E_per./(2*len)==fs(xsi);      
      len2=double(solve(fun2_strain,len));
      strain2=len2(2)*tan(xsi)/h;
      fun3 = @(locate) w*locate*strain2*E_per.*(len2(2)-locate)/len2(2);
      fun4 = @(locate) t*locate*strain2*E_per.*(len2(2)-locate)/len2(2);
      M_r=@(x)(integral(fun3,0,a1)+integral(fun4,a1,len2(2)))/1000;
      Mp=Mpc(xsi)-M_r(xsi);
      moment_load=[moment_load,Mp];
      moment_reaction=[moment_reaction,M_r(xsi)];
      moment_strap=[moment_strap,Mpc(xsi)];
      angle_matrix=[angle_matrix,xsi];
      xsi=xsi+(angle_c21-angle_c12)/n;
end
% 3rd case
xsi=angle_c22;
for k=1:n+1
      syms len
      fun3_strain = w*len*(len*tan(xsi)/h)*E_per/2==fs(xsi);      
      len3=double(solve(fun3_strain,len));
      strain3=len3(2)*tan(xsi)/h;
      fun5 = @(locate) w*locate*strain3*E_per.*(len3(2)-locate)/len3(2);
      M_r=@(x)integral(fun5,0,len3(2))/1000;
      Mp=Mpc(xsi)-M_r(xsi);
      moment_load=[moment_load,Mp];
      moment_reaction=[moment_reaction,M_r(xsi)];
      moment_strap=[moment_strap,Mpc(xsi)];
      angle_matrix=[angle_matrix,xsi];
      xsi=xsi+(0.4-angle_c22)/n;
end



end

