function Ta = temperature(Uu,Ul,Tave,phiu,phiq,param)
%
%  Revised Model version, January 2008
%
%   Integrates the wind shear over a hemisphere to get air temperature
%    Integration gives the temperature up to a constant 
%
%   Use the scalar Tave to properly set the mean global temperature

%    ------   For temperature on the staggered Q grid
Ud = (Uu-Ul);
A = -param.f0*param.a/param.R*param.dphi*cumsum((Ud(1:length(Ud)-1) + Ud(2:length(Ud)))/2);
Ta = A - sum(cos(phiq).*A)./sum(cos(phiq)) + Tave;