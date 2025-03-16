% actuation_force_4
%
% Syntax: actuation_force_4(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the actuation force 4 described in the
% paper
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function force = actuation_force_4(k,t,q,B1T,B2T,B1B,B2B,p)
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    p4 = p(4);

    force = p1*k/2*(-0.2*sin(t*2*pi*p4))*(B1T+B2T*q) + ...
            p2*k/2*(-0.2*sin(t*2*pi*3*p4))*(B1T+B2T*q) + ...
            p3*k/2*(-0.2*sin(t*2*pi*5*p4))*(B1T+B2T*q) + ...
            p1*k/2*(0.2*sin(t*2*pi*p4))*(B1B+B2B*q) + ...
            p2*k/2*(0.2*sin(t*2*pi*3*p4))*(B1B+B2B*q) + ...
            p3*k/2*(0.2*sin(t*2*pi*5*p4))*(B1B+B2B*q);
end