% actuation_signal_1
%
% Syntax: actuation_signal_1(k,t,p1)
%
% Description: function producing the actuation signal of actuation force 1
%
% Last modified: 07/02/2026, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_1(k,t,p1)
    actu = -p1*k/2*sin(t*2*pi);
end