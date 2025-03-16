% actuation_signal_4
%
% Syntax: actuation_signal_4(k,t,p)
%
% Description: function producing the actuation signal 4 described in the
% paper. 
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_4(k,t,p)
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    p4 = p(4);
    actu =  p1*k/2*(-0.2*sin(t*2*pi*p4)) + ...
        p2*k/2*(-0.2*sin(t*2*pi*3*p4)) + ...
        p3*k/2*(-0.2*sin(t*2*pi*5*p4));
end