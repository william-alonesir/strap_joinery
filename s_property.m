function [fs] = s_property(ex)
% The fitting constitutive curve of strap material
%   ex represent for extension; fs represent for strap force
%   ex = linspace(0,140); 
%   The unit for fs is kN; the unit for ex is mm
  fs=-7*10^(-8)*ex.^4 + 2*10^(-5)*ex.^3 - 0.0017*ex.^2 + 0.0989*ex; %
  syms ext_initial
end

