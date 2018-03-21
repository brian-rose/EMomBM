function [solution,Uu,Ul,Us,tau,Ta,ku,Yu,Yl,Ku,Kl,pvfu,pvfl,coalbedo,Fout,Fup,Fs,Qdot,Ko,Ks,Hocean,...
    hf,vfu,vfl,mfu,mfl,kheat,Vu,Vl,omega,omegastar,omegares,psi_e,psi_star,psi_res,Rtoa,Htotal,Hatm]...
    = final_diagnostics(Qu,Ql,Tave,Ts,S,Ai,Ad,coalbedo_o,coalbedo_i,phiu,phiq,param)

%   Efficient Model version, November 2008

[Uu,Ul,Us,tau,Ta,ice,Fs,Qdot,coalbedo,edge,ku,Yu,Yl,Ku,Kl,Ko] = current_diagnostics(Qu,Ql,Tave,Ts,S,Ai,Ad,coalbedo_o,coalbedo_i,phiu,phiq,param);
J = param.J;
kl = param.kl;
if ( ku > 0 )
   [pvfu,pvfl] = pvflux(Qu,Ql,Ku,Kl,phiu,phiq,param);
else
   pvfu = zeros(size(Uu));  pvfl = zeros(size(Ul));
end
%  Compute heat, momentum and EP fluxes
[hf,vfu,vfl,mfu,mfl,kheat] = fluxes(Uu,Ul,pvfu,pvfl,Ku,Kl,phiu,phiq,param);
%  dimensional ageostrophic mean flow (use constant f0)
Vu = -(1/param.f0) *  vfu;
Vl = -(1/param.f0) *vfl + (2*param.g / param.f0 / param.p0)*tau;
omega = -param.p0/2/param.a./cos(phiq).*diff(cos(phiu).*Vu)./diff(phiu);
%  define omegastar as the transport implied by the eddy heat flux:
omegastar = -(param.p0/2/param.Ld^2*param.R / param.f0^2 / param.a)./cos(phiq) .* diff(cos(phiu).*hf)./diff(phiu);
%  residual mean vertical motion
omegares = omega + omegastar;
%    Overturning mass streamfunctions in units of "mass Sverdrups"
psi_e = -2*pi*param.a^2 / param.g * cumtrapz(phiq,cos(phiq).*omega) * 1E-9;
psi_star = -2*pi*param.a^2 / param.g * cumtrapz(phiq,cos(phiq).*omegastar) * 1E-9;
psi_res = -2*pi*param.a^2 / param.g * cumtrapz(phiq,cos(phiq).*omegares) * 1E-9;
%  Total atmosphere + ocean heat transport implied by radiation imbalance
%     in units of W
Fup = param.Aup + param.Bup*(Ts-Ta);
Fout = param.Aout +param.Bout*Ta;
Rtoa = coalbedo.*S - Fout;
Htotal = 2*pi*param.a^2 .* [0; cumsum(cos(phiq).*Rtoa)] * param.dphi;
%  Oceanic heat transport
Tsx = [Ts(1); Ts; Ts(end)];
Ks = Ko;
if (~isempty(ice))
    Ks(ice+1) = 0;  % no diffusivity over the ice... have to add 1 because the grids are staggered!
    Tsx(min(ice+1)) = param.freeze;
end
Hocean = -2*pi.*cos(phiu)*param.Co.*Ks.*diff(Tsx)./param.dphi;
%  Atmospheric transport calculated from heating
Hatm = 2*pi*param.a^2 .* [0; cumsum(cos(phiq).*Qdot)*param.dphi];

%  This struct contains the full solution with all the diagnostics returned
%  by this function
solution.Qu = Qu; solution.Ql = Ql; solution.Tave = Tave; solution.Ts = Ts; 
solution.Uu = Uu; solution.Ul = Ul; solution.Us = Us; solution.tau = tau; solution.Ta = Ta; solution.ku = ku;
solution.Yu = Yu;  solution.Yl = Yl;  solution.Ku = Ku;  solution.Kl = Kl;
solution.pvfu = pvfu; solution.pvfl = pvfl; solution.coalbedo = coalbedo;
solution.Fout = Fout; solution.Fup = Fup; solution.Fs = Fs; solution.Qdot = Qdot;
solution.hf = hf; solution.vfu = vfu; solution.vfl = vfl; solution.mfu = mfu; solution.mfl = mfl; solution.kheat = kheat; 
solution.Vu = Vu; solution.Vl = Vl; solution.omega = omega; solution.omegastar = omegastar; solution.omegares = omegares; 
solution.psi_e = psi_e; solution.psi_star = psi_star; solution.psi_res = psi_res; solution.Rtoa = Rtoa; 
solution.Htotal = Htotal; solution.Hatm = Hatm; solution.Hocean = Hocean;
solution.S = S; solution.phiu = phiu; solution.phiq = phiq; solution.param = param; 
solution.Ko = Ko;  solution.Ks = Ks;