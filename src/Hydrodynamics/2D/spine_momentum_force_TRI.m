% spine_momentum_force_TRI3
%
% Synthax:
% tensor = spine_momentum_force_TRI(Assembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde,u,ud,udd)
%
% Description: 
% This function computes the  change of momentum along the fish spine at 
% the Assembly level (FOM), as a force (dimension is nDoFs x 1)
%
% INPUTS
%   - Assembly:             Assembly from YetAnotherFEcode.
%   - spineElementWeights:  array of length nElements, with 1 if element is
%                           a spine element, 0 else
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
%   - u:                    vector of displacements
%   - ud:                   vector of velocities
%   - udd:                  vector of accelerations
%
%
% OUTPUTS
%   - force:                change of momentum along the fish spine                  
%       
%
% Additional notes: 
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zurich
function force = spine_momentum_force_TRI(Assembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde,u,ud,udd)  
    force = Assembly.vector_hydro3('spine_momentum_force', 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, mTilde,u,ud,udd); 

end
