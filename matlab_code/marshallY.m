function Y = marshallY(U1,U3)
%
%  function Y = marshallY(U1,U3)
%
%  For the two-layer QG model, Y(y) is a function that gives the meridional
%  variation of the transfer coefficients K1 and K3.  Using the Marshall
%  (1981) formulation, it is computed from the velocities
%  U1 and U3 using the formula Y(y) = |u1-u3| / |u1-u3| max
%
%  Brian Rose
%   July 25, 2005
%
%   Note that it's a ratio, it doesn't care whether the velocities are
%   dimensional or not.
% 

Ud = abs(U1 - U3);
%dphi = pi/2/(length(U1)-1);  phiu = [0:dphi:pi/2]';
%Ud = abs(U1-U3).*cos(phiu);
m = max(Ud);
if (m > 0)
    Y = Ud / m;
else
    Y = zeros(size(Ud));
%     display('Shear is zero everywhere... cannot compute Y(y)')
end