function [solution,iceedgelat] = EMomBM(initial,S,phiu,phiq,param)
% %%%%%  EMomBM (efficient implicit scheme)
%
%     The coupled Energy Momentum Balance Model
%
%    Brian Rose
%     August 20, 2008   

J = param.J;  kl = param.kl;

if (param.numsteps>0)
     iceedgelat = zeros(1,param.numsteps);
end

%  Fix the coalbedo
coalbedo_o = param.a0 + param.a2/2*(3*sin(phiq).^2 - 1);  %ice-free, North's formula for variations with latitude
coalbedo_i = param.ai*ones(size(phiq));  % icy

%Qmodfactor = param.delt*param.R/param.Ca/param.Ld^2/param.f0./(0.5*(1-0.5*cos(phiq)));
Qmodfactor = param.delt*param.R/param.Ca/param.Ld^2/param.f0;
taufactor = 2*param.delt*param.g/param.p0/param.a./cos(phiq)./diff(phiu);
Tavefactor = param.delt/param.Ca/sum(cos(phiq));
%Kofactor = param.a*param.mu/param.Co/param.f0./cos(phiu);

%  Tridiagonal matrices Ai, Ad for the QGPV inversion, passed to fastSolveU.m
%    (only need to evaluate this once at the start, since it is fixed
%    throughout)
A1 = zeros(J+1,3);
A1(:,1) = [1; -(sec(phiq(2:J))+sec(phiq(1:J-1))).*cos(phiu(2:J)); 1];
A1(1:J-1,2) = sec(phiq(1:J-1)).*cos(phiu(1:J-1));
A1(3:J+1,3) = sec(phiq(2:J)).*cos(phiu(3:J+1));
Ai = spdiags(A1,[0,-1,1],J+1,J+1);
%const = -2/param.Ld^2*param.a^2*param.dphi^2./(0.5*(1-0.5*cos(phiu(2:end-1))));
const = -2/param.Ld^2*param.a^2*param.dphi^2;
Ad = Ai + spdiags([0; const.*ones(J-1,1); 0],0,J+1,J+1);

done = false;
count = 0;
Quold = initial.Qu; Qlold = initial.Ql; Taveold = initial.Tave; Tsold = initial.Ts;

%  begin time loop
while (~done)
    %  Here we use a timesplitting approach:  do two half-steps, first the
    %  radiation and then the diffusion
    %  wind
    [Uu,Ul,Us,tau,Ta,ice,Fs,Qdot,coalbedo,edge,ku,Yu,Yl,Ku,Kl,Ko] = current_diagnostics(Quold,Qlold,Taveold,Tsold,S,Ai,Ad,coalbedo_o,coalbedo_i,phiu,phiq,param);
    %[Uu,Ul] = fastSolveU(Quold,Qlold,Ai,Ad,phiu,phiq,param);
    %[Us,tau] = stress(Uu,Ul,param);
    % %  atmospheric temperature
    %Ta = temperature(Uu,Ul,Taveold,phiu,phiq,param);
    %[ice,Fs,Qdot] = heating(Ta,Tsold,S,coalbedo_o,coalbedo_i,param);
    %edge = min(ice);
    Tavehalf = Taveold + Tavefactor * sum(cos(phiq).*Qdot);
    Tshalf = Tsold + param.delt/param.Co*Fs;
    Quhalf = Quold - Qmodfactor.*Qdot;
    Qlhalf = Qlold + Qmodfactor.*Qdot + taufactor.*diff(cos(phiu).*tau);
    %  now do the diffusion implicitly
    Bu = Quhalf;
    Bl = Qlhalf;
    Bs = Tshalf;
    % %  PV diffusion
    %Y = marshallY(Uu,Ul);
    %Yu = Y;  Yl = Y;
    %ku = computeku(Quold,Qlold,kl,Yu,Yl,phiu);  Ku = ku*Yu;   Kl = kl*Yl;
    %if ( ku > 0 ) Ku = ku*Yu;  Kl = kl*Yl;
    %else Ku = zeros(size(phiu)); Kl = zeros(size(phiu));
    %end
    [Ku1 Ku2 Ku3] = Kvectors(Ku,phiu,phiq,param);
    [Kl1 Kl2 Kl3] = Kvectors(Kl,phiu,phiq,param);
    %%curlstress = abs(diff(cos(phiu).*tau)./diff(phiu));
    %%curlstress_interp = [0; (curlstress(1:J-1) + curlstress(2:J))/2; 0];
    %%Ko = Kofactor.*curlstress_interp;
    %%tau_interp = interp1(phiu,tau,phiq);
    %%curlstress = abs(diff(cos(phiq).*tau_interp)./diff(phiq));
    %%Ko = Kofactor.*[0;curlstress;0];
    %%Ko = param.a*param.mu/param.Co/param.f0/2/param.dphi*[0; abs( tau(3:J+1).*cos(phiu(3:J+1)) - tau(1:J-1).*cos(phiu(1:J-1)) )./cos(phiu(2:J)); 0];
    %Ko = param.a*param.mu/param.Co/param.f0/param.dphi*[0; abs(tau(2:J)-tau(1:J-1).*cos(phiu(1:J-1))./cos(phiu(2:J))); 0];
    [Ko1 Ko2 Ko3] = Kvectors(Ko,phiu,phiq,param);
    %  Atmosphere tridiagonal matrix   
    Au = sparse(diag(1+Ku2) + diag(-Ku3(1:end-1),1) + diag(-Ku1(2:end),-1));
    Al = sparse(diag(1+Kl2) + diag(-Kl3(1:end-1),1) + diag(-Kl1(2:end),-1));
    %  Ocean tridiagonal
    As = sparse(diag(1+Ko2) + diag(-Ko3(1:end-1),1) + diag(-Ko1(2:end),-1));
    %  This code modifies the diffusion problem to account for the ice edge
    if(~isempty(edge))
        %display('oh eh')
        %display(Tsold)
        if (edge==J)  % only one icy point at polar boundary
            As(edge,edge)=1;
            Bs(edge-1:edge) = Tshalf(edge-1:edge) +param.freeze*[Ko3(edge-1); -Ko1(edge)];
            As(edge-1,edge) = 0;
        elseif(edge==1)
            %  snowball...  no diffusion
            As = speye(J);
        else  % ice edge is somewhere in the interior
            As(edge+1:end,edge+1:end) = speye(length(ice)-1);
            As(edge+1,edge) = 0;  % away from ice edge we only have 1 on diagonal
            As(edge,edge:edge+1) = [1 0];
            Bs(edge-1:edge) = Tshalf(edge-1:edge) +param.freeze*[Ko3(edge-1); -Ko1(edge)];
            As(edge-1,edge) = 0;
        end
    end
    %  Time-stepping the diffusion operator is just inverting this matrix
    %  problem (sparse matrix is much more efficient at high res):
    Qunew = Au\Bu;
    Qlnew = Al\Bl;
    Tsnew = As\Bs;
    Tavenew = Tavehalf;
    
	count = count + 1;
     if (~isempty(edge))
        iceedgelat(count) = phiq(edge);
    else
        iceedgelat(count) = pi/2;
    end
    if (count==param.numsteps)
    %  Exit after a specified number of steps
        done = true;
    end
    Quold = Qunew;
    Qlold = Qlnew;
    Taveold = Tavenew;
	Tsold = Tsnew;
end
Qu = Qunew;
Ql = Qlnew;
Tave = Tavenew;
Ts = Tsnew;

%  Compute and return all diagnostics
[solution,Uu,Ul,Us,tau,Ta,ku,Yu,Yl,Ku,Kl,pvfu,pvfl,coalbedo,Fout,Fup,Fs,Qdot,Ko,Ks,Hocean,...
    hf,vfu,vfl,mfu,mfl,kheat,Vu,Vl,omega,omegastar,omegares,psi_e,psi_star,psi_res,Rtoa,Htotal,Hatm]...
    = final_diagnostics(Qu,Ql,Tave,Ts,S,Ai,Ad,coalbedo_o,coalbedo_i,phiu,phiq,param);
%[solution,Uu,Ul,Us,tau,Ta,ku,Yu,Yl,Ku,Kl,pvfu,pvfl,coalbedo,Fout,Fup,Fs,Qdot,deltaTew,Hocean,...
%    hf,vfu,vfl,mfu,mfl,kheat,Vu,Vl,omega,omegastar,omegares,psi_e,psi_star,psi_res,Rtoa,Htotal,Hatm]...
%    = final_diagnostics(Qu,Ql,Tave,Ts,S,phiu,phiq,param);