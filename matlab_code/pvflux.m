function [pvfu, pvfl] = pvflux(Qu,Ql,Ku,Kl,phiu,phiq,param);
%
%    Compute the eddy PV fluxes in each layer
%
%     For the Revised Model version, January 2008

pvgradu = diff(Qu)./diff(phiq);
pvgradl = diff(Ql)./diff(phiq);

pvfu = -Ku.*[0;pvgradu;0] / param.a;
pvfl = -Kl.*[0;pvgradl;0] / param.a;