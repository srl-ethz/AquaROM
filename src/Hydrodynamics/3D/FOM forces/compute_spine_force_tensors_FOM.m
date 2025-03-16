% compute_spine_force_tensors_FOM
%
% Synthax:
% tensor = compute_spine_force_tensors_FOM(ROMAssembly, spineElementWeights, normalisationFactors, nodeIdxPosInElements, dorsalZPos)
%
% Description: 
% This function computes the tensors that can be used to compute
% the change of momentum along the fish spine at the Assembly level (FOM).
% The tensors are then use to obtain partial derivatives used in the
% Newmark integration scheme
%
% INPUTS
%   - FOMAssembly:          Assembly from YetAnotherFEcode.
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
%   - dorsalZPos:           Position of the top of the cross-section for
%                           different spine elements. Used to compute the
%                           virtual mass of each cross-section
%
%
% OUTPUTS
%   tensors: struct with the tensors T2a, T3b, T3c, T4d (see pdf
%   documenation)
%      
% Additional notes:  -
%
% Last modified: 27/10/2024, Mathieu Dubied, ETH Zurich
function tensors = compute_spine_force_tensors_FOM(FOMAssembly, spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos)
    t0=tic;
    nel = FOMAssembly.Mesh.nElements;  % number of elements
    disp('SPINE FORCE TENSORS:')
    fprintf(' Assembling %d elements ...\n', nel)
    tic;
    
    % compute tensors
    tensors.T2a = FOMAssembly.tensor_spine_force_T2a('spine_momentum_tensor', 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
    fprintf('   2nd order term - T2a: %.2f s\n',toc)
    
    tic;
    tensors.T3b = FOMAssembly.tensor_spine_force_T3b('spine_momentum_tensor', 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos);
    tensors.T3c = FOMAssembly.tensor_spine_force_T3c('spine_momentum_tensor', 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
    fprintf('   3nd order term - T3b, T3c: %.2f s\n',toc)
    
    tic
    tensors.T4d = FOMAssembly.tensor_spine_force_T4d('spine_momentum_tensor', 'weights', spineElementWeights, nodeIdxPosInElements, normalisationFactors, dorsalZPos); 
    fprintf('   4th order term - T4d: %.2f s\n',toc)
    
    % display time needed for computation
    time = toc(t0);
    fprintf(' TOTAL TIME: %.2f s\n',time)

end
