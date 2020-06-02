function [Pc,angle_rad,k1,strain_initial] = linear_force(h,a1,a2,w,t,E_per,L)
% The first linear stage of joinery
%   Assuming no strap elongation in this stage. Pc is the critical load
%   between stage 1 and next stage
%   Unit for Pc is N, angle_rad is rad,k1 is Nm/rad
syms strain
    fun_strain = (w*a1^2+ t*(a2-a1)*(a1+a2)+w*a1*(a2+h))*strain*E_per/(2*h)==12000*0.9; %force balance
    strain_initial=double(solve(fun_strain,strain));
    fun1 = @(x) w*x*strain_initial*E_per.*(h-x)/h;
    fun2 = @(x) t*x*strain_initial*E_per.*(h-x)/h;
    Pc = (6000*0.9*h-integral(fun1,0,a1)-integral(fun2,a1,a2)-integral(fun1,a2,h))/L;   %moment balance
    angle_rad=asin(strain_initial*h/h);
    k1=Pc*L/ angle_rad/1e3;                                                             %initial stiffness
end

