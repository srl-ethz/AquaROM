% compute_actuation_tensors_FOM
%
% Synthax:
% tensors = compute_actuation_matrices_FOM(Assembly, actuationElements, rho)
%
% Description: This function computes the unreduced actuatio tensors at the
% global level, by combining, the element-level contributions. The obtained
% tensors are the ones used in the FOM.
%
% INPUTS
%  (1) Assembly:            Unreduced assembly from YetAnotherFEcode
%  (2) actuationElements:   logical array of length nElements, with 1 if an
%                           element is an actuation element, and 0 else
%  (3) actuationDirection:  contraction/extension direction
%
% OUTPUTS
%   tensors: a struct variable with the following fields:     
%       .B1                 first order tensor, i.e., a vector
%       .B2                 second order tensors, i.e., a matrix
%      	.time               computational time
%
% Last modified: 16/02/2024, Mathieu Dubied, ETH Zurich

function tensors = compute_actuation_tensors_FOM(Assembly, actuationElements, actuationDirection) 

    t0=tic;
    
    % data from ROM Assembly
    nel = Assembly.Mesh.nElements;      % number of elements
    
    % compute reduced tensors
    disp(' ACTUATION MATRICES:')
    fprintf(' Assembling %d elements ...\n', nel)

    tic;
    B1 = Assembly.vector_actuation('B1','weights', actuationElements, actuationDirection);
    fprintf('   1st order term - B1: %.2f s\n',toc)
    
    tic
    B2 = Assembly.matrix_actuation('B2','weights', actuationElements, actuationDirection);
    fprintf('   2nd order term - B2: %.2f s\n',toc)
    
    % display time needed for computation
    time = toc(t0);
    fprintf(' TOTAL TIME: %.2f s\n',time)
    
    % store outputs   
    tensors.B1 = B1;
    tensors.B2 = B2;
    tensors.time = time;

end