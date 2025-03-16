% derivatives_actuation_force_2
%
% Syntax: derivatives_actuation_force_2(k,t,q,B1T,B2T,B1B,B2B)
%
% Description: function producing the derivatives of the actuation force 2
% described in the paper
%
% Last modified: 12/05/2024, Mathieu Dubied, ETH Zürich
function der = derivatives_actuation_force_2(k,t,q,B1T,B2T,B1B,B2B,p1)
    dfdp = -k/2*0.2*cos(t*p1*pi)*pi*t*(B1T+B2T*q) + ...
            k/2*0.2*cos(t*p1*pi)*pi*t*(B1B+B2B*q);
    der.dfdp = dfdp;
end