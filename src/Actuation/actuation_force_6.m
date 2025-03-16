% actuation_force_6
%
% Syntax: actuation_force_6(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the actuation force 6 described in the
% paper
%
% Last modified: 21/01/2024, Mathieu Dubied, ETH Zürich
function force = actuation_force_6(k,t,q,B1T,B2T,B1B,B2B,p)
    p1 = p(1);  % amplitude
    p2 = p(2);  % frequency
    p3 = p(3);  % quantity of easing

    force = -(1-p3)*cos(pi*(p2*t-0.5))*(B1T+B2T*q) + ...
            - p3*cos(pi*(cos(pi*(p2*t-0.5))-1)/2)*(B1T+B2T*q) + ...
            (1-p3)*cos(pi*(p2*t-0.5))*(B1B+B2B*q) + ...
            p3*cos(pi*(cos(pi*(p2*t-0.5))-1)/2)*(B1B+B2B*q);
    force = k/2*p1*force;
end