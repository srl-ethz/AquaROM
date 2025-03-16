% all_derivatives_actuation_force_1
%
% Syntax: all_derivatives_actuation_force_1(k,t,q,B1T,B2T, B1B,B2B)
%
% Description: function producing the derivatives of the actuation force 1 
% described in the paper
%
% Last modified: 18/05/2024, Mathieu Dubied, ETH Zürich
function der = all_derivatives_actuation_force_1(k,t,q,B1T,B2T, B1B,B2B,p)
    dfdp = -k/2*sin(t*2*pi)*(B1T+B2T*q) + ...
            k/2*sin(t*2*pi)*(B1B+B2B*q);
    dfdq = -p*k/2*sin(t*2*pi)*B2T + ...
            p*k/2*sin(t*2*pi)*B2B;
    der.dfdp = dfdp;
    der.dfdq = dfdq;
end