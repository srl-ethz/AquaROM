% tail_pressure_force_TRI3
%
% Synthax:
% force = force = tail_pressure_force_TRI3(Assembly, tailElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde, q, qd)
%
% Description: 
% This function computes the tail pressure force at the Assembly level,
% using the original nonlinear formulation (not in form of a polynomial).
%
% INPUTS
%   - Assembly:             Assembly from YetAnotherFEcode.
%   - tailElementWeights:   array of length nElements, with 1 if element is
%                           the tail element, 0 else
%   - normalisationFactors: array of size nElements with zeros everywhere,
%                           except at the indexes of the spine elements.
%                           For these indexes, the inverse of the distance
%                           between the two spine nodes is entailed.
%   - nodeIdxPosInelements: matrix of size nElements x 2. For the row
%                           corresponding to spine elements, it gives the
%                           column index (1 to 3) in the element of the 
%                           node close to the tail (1st column) and the 
%                           node close to the head (2nd column).
%   - mTilde:               virtual mass linear density
%   - q:                    constrained displacement vector
%   - qd:                   derivative of the constrained displacement
%                           vector
%
% OUTPUTS
%   force:                  reactive force at the Assembly level   
%       
%
% Additional notes:
%
% Last modified: 06/10/2023, Mathieu Dubied, ETH Zurich
function force = tail_pressure_force_TRI3(Assembly, tailElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde, q, qd)
    u = Assembly.unconstrain_vector(q);
    ud = Assembly.unconstrain_vector(qd);
    force = Assembly.vector_hydro2('tail_pressure_force', 'weights', tailElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde,u,ud); 
    % force = Assembly.constrain_vector(force);
end
