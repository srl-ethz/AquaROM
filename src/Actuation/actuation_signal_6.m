% actuation_signal_6
%
% Syntax: actuation_signal_6(k,t,p)
%
% Description: function producing the actuation signal 6 described in the
% paper. 
%
% Last modified: 21/01/2024, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_6(k,t,p)
    p1 = p(1);  % amplitude
    p2 = p(2);  % frequency
    p3 = p(3);  % quantity of easing
    
    actu =  -(1-p3)*cos(pi*(p2*t-0.5)) + ...
            - p3*cos(pi*(cos(pi*(p2*t-0.5))-1)/2);
    actu = k/2*p1*actu;
end