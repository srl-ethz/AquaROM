% PROM_tail_pressure_derivates_TET4_lightweight
%
% Synthax:
% der = PROM_tail_pressure_derivatives_TET4_lightweight(q,qd,A,B,R,mTilde,w,x0,ud,VTail,UTail,z0,Uz)
%
% Description: This function returns the partial derivatives of the change
% in tail pressure force, which is understood as a force acting on the spine.
% This lightweight version only entails the partial derivatives needed in
% the algorithm proposed in the paper: dfdp for xi=0
%
% INPUTS: 
% (1) q:                vector of time domain displacements
% (2) qd:               vector of time domain velocities
% (3) A:                matrix used to express the tail pressure force
% (4) B:                matrix used to express the tail pressure force
% (5) R:                rotation matrix
% (6) w:                normalisation factor for the tail (inverse of spine
%                       section's length)

% (7) x0Tail:           initial node positions of the tail node in FOM
% (8) VTail:            part of the ROB matrix corresponding to the tail
%                       element 
% (9) UTail:            part of the shape variations matrix correponsing to
%                       the tail element
% (10) z0:              z position (undeformed) of the dorsal node matched
%                       to the tail element
% (11) Uz:              row vector containing the part of U corresponding
%                       to z0
%
% OUTPUTS:
% (1) der:              struct containing the partial derivatives
%     
%
% Additional notes:
%   - q and qd should be understood as eta and dot{eta} (ROM).
%   - p should be understood as xi.
%
% Last modified: 06/11/2023, Mathieu Dubied, ETH Zürich
function der = PROM_tail_pressure_derivatives_TET4_lightweight(q,qd,A,B,R,w,x0Tail,VTail,UTail,z0,Uz)
   
    mTilde = 0.25*pi*1000*((z0)*2)^2;
    
    % dfdp ________________________________________________________________
    firstTerm1 = A*VTail*qd;                            % vector
    firstTerm2 = dot(A*VTail*qd,R*B*(x0Tail+VTail*q));  % scalar
    firstTerm3 = B*(x0Tail+VTail*q);                    % vector
    firstTerm1RB = firstTerm1.'*R*B*UTail;              % row vector 
    
    firstTerm = 2*firstTerm2*firstTerm3*firstTerm1RB;   % outer-product between the last two terms to get a matrix
    secondTerm = dot(A*VTail*qd,R*B*(x0Tail+VTail*q))^2*B*UTail;


    fTail_wo_factor = VTail.'*(dot(A*VTail*qd,R*B*(x0Tail+VTail*q))).^2* ...
                        B*(x0Tail+VTail*q);
    thirdTerm = 0.25*pi*1000*w^3*(2*z0)*fTail_wo_factor*Uz; % specific to TET4/3D structures
    

    dfdp = 0.5*mTilde*w^3*VTail.'*(firstTerm + secondTerm) + thirdTerm;

       
    % store results in output struct ______________________________________
    der.dfdp = dfdp;      

end
