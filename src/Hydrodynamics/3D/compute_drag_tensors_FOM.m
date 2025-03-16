% compute_drag_tensors_FOM
%
% Synthax:
% tensors = compute_drag_tensors_FOM(Assembly, skinElements, skinElementFaces, rho)
%
% Description: This function computes the unreduced drag
% tensors at the global level, by combining, the
% element-level contributions. The obtained tensors are the ones used in
% the FOM.
%
% INPUTS
%   - Assembly:             Unreduced assembly from YetAnotherFEcode
%   - skinElements:         logical array of length nElements, with 1 if an
%                           element is a skin element, and 0 else
%   - skinElementFaces:     matrix of size nElements x 2, giving for each
%                           skin element which face (1,2 or 3) is part of
%                           the skin. If two faces are part of the skin, we
%                           use the 2 columns of the matrix.
%   - rho:                  water density
%
% OUTPUTS
%   tensors: a struct variable with the following fields:     
%       .T3
%      	.time           computational time
%     
%
% Additional notes:
%   - List of currently supported elements: TRI3, TET4
%
% Last modified: 18/02/2024, Mathieu Dubied, ETH Zurich

function tensors = compute_drag_tensors_FOM(Assembly, skinElements, skinElementFaces, rho) 

    t0=tic;
    
    % data from ROM Assembly
    nel = Assembly.Mesh.nElements;      % number of elements
    nDOFs = Assembly.Mesh.nDOFs;
    
    % compute reduced tensors
    disp(' DRAG TENSORS:')
    fprintf(' Assembling %d elements ...\n', nel)

    tic;
    T3 = Assembly.vector_skin('Te1','weights', skinElements, skinElementFaces, rho);

    fprintf('   3rd order term - T3: %.2f s\n',toc)
    
    % display time needed for computation
    time = toc(t0);
    fprintf(' TOTAL TIME: %.2f s\n',time)
    
    % store outputs   
    tensors.T3 = T3;
    tensors.time = time;

end


