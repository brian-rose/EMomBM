function ku = computeku(Qu,Ql,kl,Yu,Yl,phiu)
%
%   Computes the upper coefficient ku according to the integral constraint
%
%   The output k1 has same dimensions as input k3
%
%    Revised Model version, January 2008

%  Compute mean gradients
gradQubar = sum(cos(phiu).*Yu.*[0; diff(Qu); 0]);
gradQlbar = sum(cos(phiu).*Yl.*[0; diff(Ql); 0]);

if ( gradQubar ~= 0 )
    ku = - kl *  ( gradQlbar ./ gradQubar );
else
    ku = 0;
end