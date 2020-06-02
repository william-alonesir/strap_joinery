function [fs_t] = fs_t(angle)
% The relationship between the top strap force and rotation
%   The rotation angle should be within the range (0.0.6)
%   The unit for angle is rad, for ext_t is mm, for fs_t is kN
ext_t=800-200*(sqrt(10+6*cos(angle)));
fs_t=s_property(48.7-ext_t);
end