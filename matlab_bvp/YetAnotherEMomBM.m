function solution = YetAnotherEMomBM(xs,param,initial_ppy)
%   
%   Steady-state version of the Energy-Momentum Balance Model
%    Brian Rose
%   July 2 2010
%
%   This version has no ocean at all... simple EBM plus the wind stress
%   calculation
%
%   This code solves a steady-state form of the model published in Rose and Marshall, JAS 2009
%   Except the ocean component is set to zero.
%   It's an example of how to set up the steady-state boundary value
%   problem.
%    When input param.solvetau is false, the code is just solving a
%    diffusion equation for temperature, with non-linear ice-albedo feedback.
%    It is then a numerical generalization of the analytical method given in North, JAS 1975
%    "Analytical Solutions to a Simple Climate Model with Diffusive Heat Transport".
%   Coupling in the momentum budget gives another 2nd-order boundary value
%   problem for wind stress tau plus an unknown parameter.
%
%   See section 3.3 of Rose, B.E.J. (2010) "Oceanic Control of the Sea Ice
%   Edge and Multiple Equilibria in the Climate System", PhD thesis, MIT.

if nargin==3; useppyinit = true; else useppyinit = false; end

eps = 0.000001;
options = bvpset('NMax',4000);

ltheta=param.ltheta; 
A1=param.A1;  s2=param.s2;
Tf=param.Tf;  
ai=param.ai;
if (isscalar(param.ao))
    a0=param.ao;  a2=0;
elseif (length(param.ao)==2)
    a0 = param.ao(1);  a2=param.ao(2);
else
    error('input argument ao should be length 2 or scalar')
end

if (xs>=1)
    xinit = linspace(0,1-eps,50);
    modelhandle = @EMomBM_noice;
    bchandle = @EMomBMbc_noice;
    if (xs>1)
        q=param.qinitial;  % fixed q in this case
        unknownparams = false;
    else
        paramguess=param.qinitial;
        unknownparams=true;
    end
else  % a two-sided problem
    xinit = [linspace(0,xs,20) linspace(xs,1-eps,20)];
    modelhandle = @EMomBM_2sided;
    bchandle = @EMomBMbc_2sided;  
    paramguess = param.qinitial;
    unknownparams=true;
end

if (param.solvetau)
    if (xs>1)
        paramguess = param.kappaguess;
    else
        paramguess = [paramguess param.kappaguess];
    end
    unknownparams=true;
    ltau = param.ltau;
    ld = param.ld;
    s = param.s;
end

inithandle = @EMomBMinit;
if useppyinit; inithandle = @EMomBMinit_ppy; end  %override the default with user-supplied initial condition
if (unknownparams)
    solinit = bvpinit(xinit,inithandle,paramguess);
else
    solinit = bvpinit(xinit,inithandle);
end
sol = bvp4c(modelhandle,bchandle,solinit,options);

%  some post-processing
solution.sol = sol;
solution.xs = xs;
if (xs>=1)
    solution.x = linspace(0,1-eps,300);
    solution.q = q;
else
    solution.x = [linspace(0,xs-eps/2,150) linspace(xs+eps/2,1-eps,150)]; 
    solution.q = sol.parameters(1);
end
solution.lat = asin(solution.x)/pi*180;
solution.param = param;
[solution.y solution.yprime]=deval(sol,solution.x);
solution.thetastar = solution.y(1,:);
solution.Hastar = -solution.y(2,:);
if (param.solvetau)
    solution.kappa = sol.parameters(end);
    solution.taustar = solution.y(3,:)./(1-solution.x.^2);
    solution.curltaustar = solution.y(4,:);
end
for n=1:size(solution.y,1)
    solution.ppy(n) = spline(solution.x,solution.y(n,:));
end
     
    function dydx = EMomBM_noice(x,y,parameters)
        if (xs==1)
            q=parameters(1);
        end
        as = (a0+a2*P2(x))*(1+s2*P2(x));
        dsolar=q*((a0*s2+a2)*P2prime(x)+2*a2*s2*P2prime(x).*P2(x));
        dydx = [y(2,:)./(1-x.^2)
            (1/ltheta)*(y(1,:)-q*as+A1)];
        if (param.solvetau)
            kappa = parameters(end);
            dydx = [dydx; taudydx(x,y,kappa,dsolar)];
        end
    end
    function dydx = EMomBM_2sided(x,y,region,parameters)
      q = parameters(1);
      switch region
          case 1
            as = (a0+a2*P2(x))*(1+s2*P2(x));
            dsolar=q*((a0*s2+a2)*P2prime(x)+2*a2*s2*P2prime(x).*P2(x));
          case 2
            as = ai*(1+s2*P2(x));
            dsolar=q*ai*P2prime(x);
      end
      dydx = [y(2,:)./(1-x.^2)
            (1/ltheta)*(y(1,:)-q*as+A1)];
      if (param.solvetau)
            kappa = parameters(end);
            dydx = [dydx; taudydx(x,y,kappa,dsolar)];
      end
  end
    function more = taudydx(x,y,kappa,dsolar)
        more = [y(4,:)
            1+y(3,:)./ltau./(1-x.^2)+s*((1-3*kappa+(1-7*kappa)*ld/ltheta).*y(2,:)./(1-x.^2)-(1-7*kappa)*ld/ltheta*dsolar)
            y(3,:)./sqrt(1-x.^2)];
    end
    function res = EMomBMbc_2sided(YL,YR,parameters)
        q=parameters(1);
        res = [YL(2,1)  % no atm heat flux at equator
            YR(2,2)./eps  % no atm heat flux at pole
            YR(1,1)-YL(1,2)  % atm temp. continuous at ice edge
            YR(2,1)-YL(2,2)  % atm heat flux continuous at ice edge
            YR(1,1)-Tf];  %  temp at freezing threshold at ice edge
        if (param.solvetau)
            res = [res
            YL(3,1)  % zero stress at equator
            YR(3,2).*eps  % zero stress at pole
            YL(5,1)  % zero torque integral at equator
            YR(5,2)  % zero torque integral at pole
            YR(3,1)-YL(3,2)  % stress continuous at ice edge
            YR(4,1)-YL(4,2)  % curl(stress) continuous at ice edge
            YR(5,1)-YL(5,2)];  % torque integral continuous at ice edge
        end
    end
    function res = EMomBMbc_noice(ya,yb,parameters)
        if (xs<=1)
            q=parameters(1);
        end
        res = [ya(2)  % no atm heat flux at equator
            yb(2)./eps];  % no atm heat flux at pole
        if (xs==1)
            res = [res;yb(1)-Tf];  %  temp at freezing threshold at pole
        end
        if (param.solvetau)
            res = [res
                ya(3)  % zero stress at equator
                yb(3).*eps  % zero stress at pole
                ya(5)  % zero torque integral at equator
                yb(5)];  % zero torque integral at pole
        end
    end
    function yinit = EMomBMinit_ppy(x,region,parameters)
        yinit = [ppval(initial_ppy(1),x)
            ppval(initial_ppy(2),x)];
        if (param.solvetau)
            yinit = [yinit
                ppval(initial_ppy(3),x)
                ppval(initial_ppy(4),x)
                ppval(initial_ppy(5),x)];
        end
    end
    function yinit = EMomBMinit(x,region,parameters)
        yinit = [0.5-P2(x)
            -(1-x.^2).*x*3];
        if (param.solvetau)
            yinit = [yinit
                -sin(x)
                cos(x)
                -cos(x)];
        end
    end
        function y = P2(x)
        %  the second order Legendre polynomial
            y = 1/2*(3*x.^2-1);
        end
        function y = P2prime(x)
            y = 3*x;
        end
    function y = P4(x)
        %  the fourth order Legendre polynomial
            y = 1/8*(35*x.^4-30*x.^2+3);
    end
    function S = solar(x)
    %  This is a very accurate fit to the annual mean aquaplanet solar
    %  forcing (divided by S0/4)
		S = (1.2371-0.6810*x.^2-0.3490*x.^4+1.2320*x.^6-1.7483*x.^8+0.8579*x.^10);
    end
end