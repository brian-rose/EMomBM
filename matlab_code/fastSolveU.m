function [U1,U3] = fastSolveU(Q1,Q3,Ai,Ad,phiu,phiq,param)
%
%  function [U1,U3] = fastSolveU(Q1,Q3,Ai,Ad,phiu,phiq,param)
%
%     For the two-layer quasi-geostrophic model, solves for the zonal
%     velocities U from the potential vorticities Q
%  Layer 1 is the upper layer, Layer 3 is the lower layer.
%
%   Updated for the Efficient Model, August 2008 --- using sparse matrix
%   inversion for speed, and no need to recalculate the LHS tridiagonal
%   matrix for each iteration...  simply compute those at beginning of time
%   loop and pass to this function.

%  Compute sum and difference:
Qd = Q1 - Q3;
Qi = Q1 + Q3;

%  This code creates the RHS of the equations for Uitilda and Udtilda
%Bi = [0; param.a*param.dphi*(4*param.omega*diff(sin(phiq)) - diff(Qi)); 0];
Bi = [0; param.a*param.dphi*(4*param.f0/sqrt(2)*diff(sin(phiq)) - diff(Qi)); 0];
Bd = -param.a*param.dphi*[0; diff(Qd); 0];

%  Solve for Ui (the depth integral) and Ud (the difference)
Uitilda = Ai \ Bi;
Udtilda = Ad \ Bd;

% %  Get Ui and Ud by dividing out the factor cos(latitude)
%Ui = [0; Uitilda(2:param.J) ./ cos(phiu(2:param.J)); 0];
%Ud = [0; Udtilda(2:param.J) ./ cos(phiu(2:param.J)); 0];
Ui = [0; Uitilda(2:param.J); 0];
Ud = [0; Udtilda(2:param.J); 0];

%  Solve for U1 and U3
U1 = (1/2)*(Ui+Ud);
U3 = (1/2)*(Ui-Ud);