% compute_spine_momentum_tensor_PROM_TET4
%
% Synthax:
% tensor = compute_spine_momentum_tensor_PROM_TET4(PROMAssembly, spineElementWeights, normalisationFactors, nodeIdxPosInElements, mTilde)
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
% Last modified: 22/10/2023, Mathieu Dubied, ETH Zurich
function tensors = compute_spine_momentum_tensor_PROM_TET4(PROMAssembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx)
    m = size(PROMAssembly.V,2);
    md = size(PROMAssembly.U,2);

    % f0, f1, f2 (tensors.Txx.f0, tensors.Txx.f1, tensors.Txx.f2)
    tensors.Txx = PROMAssembly.tensor_spine_momentum_xx_TET4_PROM('spine_momentum_tensor',[m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx); 
    tensors.TxV = PROMAssembly.tensor_spine_momentum_xV_TET4_PROM('spine_momentum_tensor',[m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx); 
    tensors.TVx = PROMAssembly.tensor_spine_momentum_Vx_TET4_PROM('spine_momentum_tensor',[m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx); 
    tensors.TVV = PROMAssembly.tensor_spine_momentum_VV_TET4_PROM('spine_momentum_tensor',[m m m m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx); 
    tensors.TUV = PROMAssembly.tensor_spine_momentum_UV_TET4_PROM('spine_momentum_tensor',[m m md m], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx); 
    tensors.TVU = PROMAssembly.tensor_spine_momentum_VU_TET4_PROM('spine_momentum_tensor',[m m m md], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx); 
    tensors.TUU = PROMAssembly.tensor_spine_momentum_UU_TET4_PROM('spine_momentum_tensor',[m m md md], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx);
    tensors.TxU = PROMAssembly.tensor_spine_momentum_xU_TET4_PROM('spine_momentum_tensor',[m m md], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx);
    tensors.TUx = PROMAssembly.tensor_spine_momentum_Ux_TET4_PROM('spine_momentum_tensor',[m m md], 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos, dorsalNodeIdx);

end
