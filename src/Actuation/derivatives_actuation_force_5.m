% derivatives_actuation_force_5
%
% Syntax: derivatives_actuation_force_5(k,t,q,B1T,B2T,B1B,B2B,p)
%
% Description: function producing the derivatives of the actuation force 5
% described in the paper
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zürich
function der = derivatives_actuation_force_5(k,t,q,B1T,B2T,B1B,B2B,p)
    h = 0.005;
    timeStep = int32(t/h)+1;
    gradientAtTimeStep = zeros(1,length(p));
    gradientAtTimeStep(timeStep) = 1;

    dfdp = (k/2*0.2*(B1T+B2T*q) + k/2*-0.2*(B1B+B2B*q))*gradientAtTimeStep;
    der.dfdp = dfdp;
end