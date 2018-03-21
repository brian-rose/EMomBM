function [hf,vfu,vfl,mfu,mfl,kheat] = fluxes(Uu,Ul,pvfu,pvfl,Ku,Kl,phiu,phiq,param)
%
%   computes the heat and relative vorticity fluxes
%  
%    Revised Model version, January 2008

Y = Kl/param.kl;
%  magnitude of heat flux parameter, in accordance with integral constraint
kheat = -param.Ld^2 * trapz(phiu,cos(phiu).*pvfu) ./ trapz(phiu,cos(phiu).*Y.*(Uu - Ul));

%  heat flux
hf = (param.f0 / param.R) * kheat * Y.*(Uu-Ul);

%  relative vorticity fluxes
vfu = pvfu + param.R/param.Ld^2/param.f0*hf;
vfl = pvfl - param.R/param.Ld^2/param.f0*hf;

%  momentum fluxes...  integrate across hemisphere
mfu = -param.a./cos(phiu).*cumtrapz(phiu, cos(phiu).*vfu); 
mfl = -param.a./cos(phiu).*cumtrapz(phiu, cos(phiu).*vfl); 
%  force them to zero at the pole
mfu(param.J+1) = 0;
mfl(param.J+1) = 0;