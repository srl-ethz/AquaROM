% derivatives_actuation_force_3
%
% Syntax: derivatives_actuation_force_3(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the derivatives of the actuation force 3
% described in the paper
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function der = derivatives_actuation_force_3(k,t,q,B1T,B2T,B1B,B2B,p)
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);

    dfdp = [k/2*(-0.2*sin(t*p2+p3))*(B1T+B2T*q) + ...
            k/2*(0.2*sin(t*p2+p3))*(B1B+B2B*q), ...
            p1*k/2*(-0.2*cos(t*p2+p3))*t*(B1T+B2T*q) + ...
            p1*k/2*(0.2*cos(t*p2+p3))*t*(B1B+B2B*q), ...
            p1*k/2*(-0.2*cos(t*p2+p3))*(B1T+B2T*q) + ...
            p1*k/2*(0.2*cos(t*p2+p3))*(B1B+B2B*q)];
    der.dfdp = dfdp;
end