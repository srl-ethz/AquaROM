% actuation_force_1
%
% Syntax: actuation_force_1(k,t,q,B1T,B2T,B1B,B2B,p1)
%
% Description: function producing the actuation force 1 described in the
% paper
%
% Last modified: 12/05/2024, Mathieu Dubied, ETH Zürich
function force = actuation_force_1(k,t,q,B1T,B2T,B1B,B2B,p1)
    force = -p1*k/2*sin(t*2*pi)*(B1T+B2T*q) + ...
            p1*k/2*sin(t*2*pi)*(B1B+B2B*q);
end