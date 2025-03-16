% tensors_hydro_FOM
%
% Synthax:
% tensors = tensors_hydro_FOM(myAssembly, elements, skinElements, skinElementFaces, vwater, rho)
%
% Description: 
% This function computes the hydrodynamic tensors at the
% Assembly level, by combining the element-level tensors. It uses classical
% element-to-assembly combining method, and results in FOM (large) tensors.
%
%
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - elements: table of the elements
%   - skinElements:
%   - skinElementFaces
%   - vwater
%   - rho:
%   - c: scaling factor for thrust force
%
% OUTPUT:
%   tensors: a struct variable with the following fields*:
%       .T1
%       .Tu2
%       .Tudot2
%       .Tudot3
%       .Tuu3
%       .Tuudot3
%    	.Tudotudot3         
%      	.time           computational time
%     
%
% Additional notes:
%   - ALL the elements are assumed to have the same properties in terms
%     of MATERIAL and QUADRATURE rules.
%   - List of currently supported elements: 
%     TRI3
%
% Last modified: 16/04/2023, Mathieu Dubied, ETH Zurich

function tensors = unreduced_tensors_hydro_FOM(myAssembly, elements, skinElements, skinElementFaces, vwater, rho,c)

t0=tic;
% data from myAssembly
%nodes    = myAssembly.Mesh.nodes;                   % nodes table
nel      = myAssembly.Mesh.nElements;           	% number of elements
%nnodes   = myAssembly.Mesh.nNodes;               	% number of nodes
%freedofs = myAssembly.Mesh.EBC.unconstrainedDOFs;   % free DOFs
nDOFs = myAssembly.Mesh.nDOFs;

myMesh = myAssembly.Mesh;
u0 = zeros(myMesh.nDOFs, 1);    


% compute FOM tensors
disp(' FOM HYDRODYNAMIC TENSORS:')
fprintf(' Assembling %d elements ...\n', nel)

tic
T1 = myAssembly.vector_skin('Te1', 'weights', skinElements, skinElementFaces, vwater, rho,c);
fprintf('   1st order terms - T1: %.2f s\n',toc)

tic
Tu2 = myAssembly.matrix_skin('Teu2', 'weights', skinElements, skinElementFaces, vwater, rho,c);
Tudot2 = myAssembly.matrix_skin('Teudot2', 'weights', skinElements, skinElementFaces, vwater, rho,c);
fprintf('   2nd order terms - Tu2, Tudot2: %.2f s\n',toc)


tic;
Tuu3 = 0.5*myAssembly.tensor_skin('Teuu3',[nDOFs nDOFs nDOFs],[2 3],'weights', skinElements, skinElementFaces, vwater, rho,c);
Tuudot3 = myAssembly.tensor_skin('Teuudot3',[nDOFs nDOFs nDOFs],[2 3], 'weights', skinElements, skinElementFaces, vwater, rho,c);
Tudotudot3 = 0.5*myAssembly.tensor_skin('Teudotudot3',[nDOFs nDOFs nDOFs],[2 3], 'weights', skinElements, skinElementFaces, vwater, rho,c);
fprintf('   3rd order terms - Truu3, Truudot3, Turudotudot3: %.2f s\n',toc)


time = toc(t0);

fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)


% outputs
tensors.T1 = T1;
tensors.Tu2 = Tu2;
tensors.Tudot2 = Tudot2;
tensors.Tuu3 = Tuu3;
tensors.Tuudot3 = Tuudot3;
tensors.Tudotudot3 = Tudotudot3; 
tensors.time = time;

end


