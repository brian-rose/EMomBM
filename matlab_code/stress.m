function [Us, tau] = stress(Uu,Ul,param)
%
%   Revised Model version, January 2008
%   Computes the wind at the surface Us and the stress ON the surface tau
%   (the wind stress)

Us = (3/2)*Ul - (1/2)*Uu;
tau = param.eps * Us;