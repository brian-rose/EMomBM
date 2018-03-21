function [ice_indices,Fs,Qdot,coalbedo] = heating(Ta,Ts,S,coalbedo_o,coalbedo_i,param)
%  Calculate the atmosphere-ocean heat exchange terms
%  This function is called by the full model as well as the 
%  Two-Level EBM without dynamics
%
%    For the Effifient Model, August 2008
%
%    This version is designed for a single surface temperature Ts...
%    represents ice surface where cold, and ocean surface otherwise

coalbedo = coalbedo_o;  %ice-free
%  Albedo depends on ice state:
ice_indices = find(Ts<param.freeze);
if (~ isempty(ice_indices) )
    coalbedo(ice_indices) = coalbedo_i(ice_indices);
end
%  Outgoing longwave
Fout = param.Aout + param.Bout*Ta;
%  Net surface to atmosphere flux
Fup = param.Aup + param.Bup*(Ts-Ta);
%  this is an experiment
Fup(ice_indices) = param.Aup + (param.Bup-7)*(Ts(ice_indices)-Ta(ice_indices));
Fs = coalbedo.*S - Fup;  % net heat flux IN to surface in W/m^2
Qdot = Fup-Fout;   % net atmospheric heating in W/m^2