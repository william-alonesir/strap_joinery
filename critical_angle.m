function [angle_c11,angle_c21,angle_c12,angle_c22] = critical_angle(angle_rad,h,a1,a2,w,t,E_per,L)
% Obtain two critical strain in this stage
%   The strain_c1 is the critical value between three area affected and two
%   area affected. The strain_c2 is the critical value between two areas
%   affected and one area affected

syms lenc1 lenc2
n1=0;n2=0;
alpha1=@(x)atan((L*sin(x)-h*cos(x)+h)/(L*cos(x)+h*sin(x)+h));       % in rad
alpha2=@(x)atan((L*sin(x)/(L*cos(x)+h)));                           % in rad
fs=@(x)1000*(2*fs_b(x)*cos(alpha1(x))+2*fs_t(x)*cos(alpha2(x)));    % in N
for xsi=angle_rad:0.00005:0.01
    n1=n1+1;
    
    fun_len1 = (w*(lenc1-a2).^2+t*(a2-a1)*(2*lenc1-a1-a2)+w*a1*(2*lenc1-a1)).*(tan(xsi)*lenc1/h)*E_per./(2*lenc1)==fs(xsi);
    lenc1_increase=double(solve(fun_len1,lenc1))-a2;
    if lenc1_increase(2)<0
       angle_c12=xsi;
       break    
    end
    angle_c11=xsi;
end

for xsi=angle_c12:0.001:0.6
    n2=n2+1;
    fun_len2=(t*(lenc2-a1)^2+w*a1*(2*lenc2-a1)).*(tan(xsi)*lenc2/h)*E_per./(2*lenc2)==fs(xsi);
    lenc2_increase=double(solve(fun_len2,lenc2))-a1;
    if lenc2_increase(2)<0
        angle_c22=xsi;
        break    
    end
    angle_c21=xsi;
end


end

