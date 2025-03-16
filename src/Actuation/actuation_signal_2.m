% actuation_signal_2
%
% Syntax: actuation_signal_2(k,t,p)
%
% Description: function producing the actuation signal 2 described in the
% paper. The actuation signal is defined as df_actu/dq
%
% Last modified: 12/05/2024, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_2(k,t,p1)
    actu = -k/2*0.2*sin(t*p1*pi);
end