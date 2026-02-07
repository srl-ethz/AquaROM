% actuation_force_1
%
% Syntax: actuation_force_1(k,t,q,B1T,B2T,B1B,B2B,p1)
%
% Description: function producing a simple acutation force scaled by p1, an
% optimization parameter
%
% Last modified: 07/02/2026, Mathieu Dubied, ETH Zürich
function force = actuation_force_1(k,t,q,B1T,B2T,B1B,B2B,p1)
    force = -p1*k/2*sin(t*2*pi)*(B1T+B2T*q) + ...
            p1*k/2*sin(t*2*pi)*(B1B+B2B*q);
end