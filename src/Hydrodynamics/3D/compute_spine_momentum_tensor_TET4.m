% compute_spine_momentum_tensor_TET4
%
% Synthax:
% tensor = compute_spine_momentum_tensor_TET4(ROMAssembly, spineElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde)
%
% Description: 
% This function computes the 4th order tensor that can be used to compute
% the change of momentum along the fish spine at the Assembly level (FOM)
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
%
%
% OUTPUTS
%   tensor: the 4th order tensor mentioned above
%       
%
% Additional notes:  
%
% Last modified: 19/10/2023, Mathieu Dubied, ETH Zurich
function tensors = compute_spine_momentum_tensor_TET4(ROMAssembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos)
    m = size(ROMAssembly.V,2);
    tensors.Txx = ROMAssembly.tensor_spine_momentum_xx_TET4('spine_momentum_tensor',[m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
    tensors.TxV = ROMAssembly.tensor_spine_momentum_xV_TET4('spine_momentum_tensor',[m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
    tensors.TVx = ROMAssembly.tensor_spine_momentum_Vx_TET4('spine_momentum_tensor',[m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
    tensors.TVV = ROMAssembly.tensor_spine_momentum_VV_TET4('spine_momentum_tensor',[m m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
end
