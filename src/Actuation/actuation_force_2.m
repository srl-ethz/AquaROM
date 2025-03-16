% actuation_force_2
%
% Syntax: actuation_force_2(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the actuation force 2 described in the
% paper
%
% Last modified: 12/05/2024, Mathieu Dubied, ETH Zürich
function force = actuation_force_2(k,t,q,B1T,B2T,B1B,B2B,p1)
    force = -k/2*0.2*sin(t*p1*pi)*(B1T+B2T*q) + ...
            k/2*0.2*sin(t*p1*pi)*(B1B+B2B*q);
end