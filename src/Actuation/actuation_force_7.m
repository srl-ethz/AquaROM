% actuation_force_7
%
% Syntax: actuation_force_7(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the actuation force 7 described in the
% paper
%
% Last modified: 12/05/2024, Mathieu Dubied, ETH Zürich
function force = actuation_force_7(k,t,q,B1T,B2T,B1B,B2B,p)
    p1 = p(1);  % amplitude
    p2 = p(2);  % frequency
   
    force = -sin(pi*(2+p2)*t)*(B1T+B2T*q) + ...
           sin(pi*(2+p2)*t)*(B1B+B2B*q);
    force = k/2*(p1+0.2)*force;
end