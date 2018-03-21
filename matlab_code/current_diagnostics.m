function    [Uu,Ul,Us,tau,Ta,ice,Fs,Qdot,coalbedo,edge,ku,Yu,Yl,Ku,Kl,Ko] = current_diagnostics(Qu,Ql,Tave,Ts,S,Ai,Ad,coalbedo_o,coalbedo_i,phiu,phiq,param)
%
%   Brian Rose
%    Efficient Model version, August 2008
%
J = param.J;
[Uu,Ul] = fastSolveU(Qu,Ql,Ai,Ad,phiu,phiq,param);
[Us,tau] = stress(Uu,Ul,param);
%curlstress = abs(diff(cos(phiu).*tau)./diff(phiu));
%Ko = param.a*param.mu/param.Co/param.f0./cos(phiu).*[0; (curlstress(1:J-1) + curlstress(2:J))/2; 0];
%Ko = param.a*param.mu/param.Co/param.f0/param.dphi*[0; abs(tau(2:J)-tau(1:J-1).*cos(phiu(1:J-1))./cos(phiu(2:J))); 0];
curltau = [0; (tau(3:J+1).*cos(phiu(3:J+1)) - tau(1:J-1).*cos(phiu(1:J-1)))./cos(phiu(2:J))/param.dphi/2; 0];
Ko = param.mu*param.a/param.Co/param.f0.*curltau.^2.*cos(phiu);
Ta = temperature(Uu,Ul,Tave,phiu,phiq,param);
[ice,Fs,Qdot,coalbedo] = heating(Ta,Ts,S,coalbedo_o,coalbedo_i,param);
edge = min(ice);
Y = marshallY(Uu,Ul);
Yu = Y;  Yl = Y;
kl = param.kl;
ku = computeku(Qu,Ql,kl,Yu,Yl,phiu);
if ( ku > 0 ) Ku = ku*Yu;  Kl = kl*Yl;
else Ku = zeros(size(phiu)); Kl = zeros(size(phiu));
end

% Tsx = [Ts(1); Ts; Ts(length(Ts))];
% deltaTew = -diff(Tsx)./param.dphi.*cos(phiu)*2*pi*param.mu;   %East-west temperature contrast across ocean basin...
% %  Note the extra factor of cos(phi) in this version, which was absent in earlier model versions.
% J = param.J;
% curlstress = abs(diff(cos(phiu).*tau)./diff(phiu));
% curlstress_interp = [0; (curlstress(1:J-1) + curlstress(2:J))/2; 0];
% %   tau_interp = interp1(phiu,tau,phiq);
% %    curlstress = abs(diff(cos(phiq).*tau_interp)./diff(phiq));
% %    Kofactor = pi*param.a*param.mu*param.co/param.omega/param.Co./cos(phiu);
% %    Ko = Kofactor.*[0;curlstress;0];
% %CKo = pi*param.a*param.mu*param.co/param.omega./cos(phiu).*curlstress_interp;
% CKo = param.a*param.mu/param.f0./cos(phiu).*curlstress_interp;
% 
% %CKo = param.a*param.mu/param.f0/param.dphi*[0; abs(tau(2:J).*cos(phiu(2:J))-tau(1:J-1).*cos(phiu(1:J-1)))./cos(phiq(1:J-1)); 0];
% 
% %CKs = Ko*param.Co;
% CKs = CKo;
% CKs(ice_indices+1) = 0;  % no diffusivity over the ice... have to add 1 because the grids are staggered!
% if (~isempty(ice_indices))
%     Tsx(min(ice_indices+1)) = param.freeze;
% end
% Hocean = -2*pi.*cos(phiu).*CKs.*diff(Tsx)./param.dphi;