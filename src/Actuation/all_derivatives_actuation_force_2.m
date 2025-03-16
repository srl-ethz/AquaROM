% all_derivatives_actuation_force_2
%
% Syntax: all_derivatives_actuation_force_2(k,t,q,B1T,B2T, B1B,B2B)
%
% Description: function producing the derivatives of the actuation force 1 
% described in the paper
%
% Last modified: 18/05/2024, Mathieu Dubied, ETH Zürich
function der = all_derivatives_actuation_force_2(k,t,q,B1T,B2T, B1B,B2B,p)
    dfdp = -k/2*0.2*cos(t*p*pi)*pi*t*(B1T+B2T*q) + ...
            k/2*0.2*cos(t*p*pi)*pi*t*(B1B+B2B*q);
    dfdq = -0.2*k/2*sin(t*p*pi)*B2T + ...
            0.2*k/2*sin(t*p*pi)*B2B;
    der.dfdp = dfdp;
    der.dfdq = dfdq;
end