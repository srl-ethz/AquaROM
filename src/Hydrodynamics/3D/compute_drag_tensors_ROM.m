% compute_drag_tensors_ROM
%
% Synthax:
% tensors = compute_drag_tensors_ROM(ROMAssembly, skinElements, skinElementFaces, rho)
%
% Description: This function computes the reduced order drag
% tensors at the global level, by combining (i.e., summing), the
% element-level contributions. The obtained tensors are the ones used in
% the ROM.
%
% INPUTS
%   - ROMAssembly:          Reduced assembly from YetAnotherFEcode
%   - skinElements:         logical array of length nElements, with 1 if an
%                           element is a skin element, and 0 else
%   - skinElementFaces:     matrix of size nElements x 2, giving for each
%                           skin element which face (1,2 or 3) is part of
%                           the skin. If two faces are part of the skin, we
%                           use the 2 columns of the matrix.
%   - rho:                  water density
%   - VHead:                row of V corresponding to the head x DOF
%
% OUTPUTS
%   tensors: a struct variable with the following fields:     
%       .Tr3
%      	.time           computational time
%     
%
% Additional notes:
%   - List of currently supported elements: TRI3, TET4
%
% Last modified: 27/10/2023, Mathieu Dubied, ETH Zurich

function tensors = compute_drag_tensors_ROM(ROMAssembly, skinElements, skinElementFaces, rho, VHead) 

    t0=tic;
    
    % data from ROM Assembly
    nel = ROMAssembly.Mesh.nElements;      % number of elements
    V = ROMAssembly.V;
    m = size(V,2);                          % size of the ROM
    
    % compute reduced tensors
    disp(' REDUCED HYDRODYNAMIC TENSORS:')
    fprintf(' Assembling %d elements ...\n', nel)

    tic;
    Tr3 = ROMAssembly.tensor_skin('Te1',[m m m],'weights', skinElements, skinElementFaces, rho, VHead);

    fprintf('   3rd order term - Tr3: %.2f s\n',toc)
    
    % display time needed for computation
    time = toc(t0);
    fprintf(' TOTAL TIME: %.2f s\n',time)
    
    % store outputs   
    tensors.Tr3 = Tr3;
    tensors.time = time;

end


