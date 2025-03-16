% actuation_signal_3
%
% Syntax: actuation_signal_3(k,t,p)
%
% Description: function producing the actuation signal 3 described in the
% paper. 
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_3(k,t,p)
    p1 = p(1);
    p2 = p(2);
    p3 = p(3);
    actu =  p1*k/2*(-0.2*sin(t*p2+p3));
end