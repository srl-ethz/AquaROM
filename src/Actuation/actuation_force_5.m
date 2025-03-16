% actuation_force_5
%
% Syntax: actuation_force_5(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the actuation force 5 described in the
% paper
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function force = actuation_force_5(k,t,q,B1T,B2T,B1B,B2B,p)
    h = 0.005;
    timeStep = int32(t/h)+1;
    
    force = p(timeStep)*k/2*0.2*(B1T+B2T*q) + ...
            p(timeStep)*k/2*-0.2*(B1B+B2B*q);
end