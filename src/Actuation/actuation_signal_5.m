% actuation_signal_5
%
% Syntax: actuation_signal_5(k,t,p)
%
% Description: function producing the actuation signal 5 described in the
% paper. 
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function actu = actuation_signal_5(k,t,p)
    h = 0.005;
    timeStep = int32(t/h)+1;
    
    actu =  p(timeStep)*k/2*0.2;
end