function [fs_b] = fs_b(angle)
% The relationship between the bottom strap force and rotation
%   The rotation angle should be within the range (0.0.6)
%   The unit for angle is rad, for ext_b is mm, for fs_b is kN
ext_b=400*(sqrt(3+cos(angle)+2*sin(angle))-2);
fs_b=s_property(48.7+ext_b);
end

