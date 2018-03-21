%%%  Script to do an iceline experiment for the EMomBM

%  Brian Rose
%  August 2008

initialS0 = 1420;
deltaS0 = 1;

param.J = 90;  % number of grid points
param.delt = 3600*24*3;  %  three day timestep
param.numsteps = 6*365*24*3600/param.delt;   % Number of timesteps ( 6 years)
param.p0 = 1E5;  % Surface pressure in units of Pa
param.f0 = 1E-4;   %  constant Coriolis parameter, in units of s^-1
param.omega = 2*pi/(24*3600);  % earth's rotation rate in rad/s
param.cp = 1004;  %  specific heat of dry air at constant pressure, in J/kg/K
param.R = 287;   %  gas constant for dry air, in J/kg/K
param.co = 4187;  % specific heat of sea water, in J / kg / K
param.g = 9.8;  %  gravity in m/s^2
param.a =  6.373E6;   %  Earth's radius in m
param.eps = 4E-2;  %  bottom stress parameter, in units of kg / m^2 / s
param.Ld = 4.4E5;   %  Rossby deformation radius for atmosphere, in units of m
param.kl = 12E6;  % magnitude of PV diffusion, in units of m^2 / s
param.s2 = -0.48;   %  solar expansion coefficient
param.freeze = -10;   %  Sea surface freezing threshhold, in degrees C
param.ai = 0.38;   %  ice coalbedo
param.a0 = 0.697;  % ice-free coalbedo expansion coefficients
param.a2 = -0.0779;
param.Aout = 220;   %  in W/m^2
param.Bout = 1.71;  %  in W/m^2/degreeC
param.Aup = 238.38;
param.Bup = 15.26;   % Upward surface heat flux parameter, in W/m^2/degreeC
param.mu = 0.012;  % Cross-basin temperature difference relative to meridional temperature difference (dimensionless)
param.dphi = pi/2 / param.J;   %  grid spacing (phi is latitude ) in units of radians
param.Ca = param.cp*param.p0/param.g;   %  column-wise heat capacity, atmosphere (J/m^2/K)
param.Co = param.Ca;

phiu = [0:param.dphi:pi/2]';
phiq = [param.dphi/2:param.dphi:pi/2-param.dphi/2]';

initial = noflow(phiu,phiq,param);

done = false;
s=1;
param.S0 = initialS0;
display(['Cooling down from initial S0 = ' num2str(initialS0)])
while( ~done)
    S0array(s) = param.S0;
    display(param.S0)
	S = param.S0/4*(1 + param.s2/2*(3*sin(phiq).^2-1));
	[solution,iceedgelat] = EMomBM(initial,S,phiu,phiq,param);
    iceline_solution(s) = solution;
    icelat(s) = iceedgelat(end);
	if(icelat(s) == min(phiq))  %  we're frozen over.. go back to the previous iteration and start warming up
        done = true;  
        display(['Found snowball at S0 = ' num2str(param.S0)])
        param.S0 = param.S0 + 2*deltaS0;
    else
        initial = iceline_solution(s);
        param.S0 = param.S0 - deltaS0;
        s = s+1;
    end
end
%  now we warm
display(['Warming up from S0 = ' num2str(param.S0)])
done = false;
while (~done)
    S0array(s) = param.S0;
    display(param.S0)
    S = param.S0/4*(1 + param.s2/2*(3*sin(phiq).^2-1));
	[solution,iceedgelat] = EMomBM(initial,S,phiu,phiq,param);
    iceline_solution(s) = solution;
    icelat(s) = iceedgelat(end);
    if (param.S0 == initialS0)
        done=true;
    else
        initial = iceline_solution(s);
        param.S0 = param.S0 + deltaS0;
        s = s+1;
    end
end

save iceline_qg_smallBout.mat 