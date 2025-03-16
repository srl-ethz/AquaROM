% actuation_signal_7
%
% Syntax: actuation_signal_7(k,t,p)
%
% Description: function producing the actuation signal 7 described in the
% paper. 
%
% Last modified: 12/05/2024, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_7(k,t,p)
    p1 = p(1);  % amplitude
    p2 = p(2);  % frequency
    
    actu =  -sin(pi*(2+p2)*t);
    actu = k/2*(0.2+p1)*actu;
end