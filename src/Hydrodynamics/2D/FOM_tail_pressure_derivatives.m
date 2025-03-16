% FOM_tail_pressure_derivates
%
% Synthax:
% der = FOM_tail_pressure_derivatives(q,qd,A,B)
%
% Description: This function returns the partial derivatives of the tail
% pressure force needed to solve the sensitity ODE as well as the EoMs.
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qd:               vector of time domain velocities
% (3) A:                matrix used to express the tail pressure force
% (4) B:                matrix used to express the tail pressure force
% (5) x0:               initial node positions of the tail node
% (6) ud:               deformed initial shape, based on shape variation
%
% OUTPUTS:
% (1) der:              strucure array containing first (and possibly
%                       second) order partial derivatives needed to solve 
%                       the sensitivity ODE and EoMs. 
%     
%
% Additional notes: -
%
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zürich
function der = FOM_tail_pressure_derivatives(q,qd,A,B,w,x0,ud)

    % dfdq ________________________________________________________________
    firstTerm1 = A*qd;                      % vector 6x1
    firstTerm2 = dot(A*qd,B*(x0+ud+q));     % scalar
    firstTerm3 = w*B*(x0+ud+q);             % vector 6x1
    firstTerm1B = (B.'*firstTerm1).';       % row vector 1x6
    
    firstTerm = 2*firstTerm2*firstTerm1B*firstTerm3;    % outer-product between the last two terms to get a matrix 6x6
    secondTerm = dot(A*qd,B*(x0+ud+q))^2*w*B;
    
    dfdq = firstTerm + secondTerm;
    
    % dfdqd _______________________________________________________________
    termB = B*(x0+ud+q);                    % vector 6x1
    term1B = (A.'*termB).';                 % row vector 1x6
    term2 = dot(A*qd,B*(x0+ud+q));          % scalar
    term3 = w*B*(x0+ud+q);                  % vector 1x6
    dfdqd = term2*term3*term1B;             % with outer-product, 6x6
    
    % store results in output struct ______________________________________
    der.dfdq = dfdq;
    der.dfdqd = dfdqd;

end
