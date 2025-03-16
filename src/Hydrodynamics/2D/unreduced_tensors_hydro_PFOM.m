% unredecued_tensors_hydro_PFOM
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
%
% OUTPUT:
%   tensors: a struct variable with the following fields*:
%       .T1
%       .T2
%       .Tu2
%       .Tu3
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
%   - no 4th order tensors for now
%
% Last modified: 16/03/2023, Mathieu Dubied, ETH Zurich

function tensors = unreduced_tensors_hydro_PFOM(myAssembly, elements, U, skinElements, skinElementFaces, vwater, rho)

t0=tic;
% data from myAssembly
nel      = myAssembly.Mesh.nElements;           	% number of elements
nDOFs = myAssembly.Mesh.nDOFs;
md = size(U,2);                         % size of the defect basis 


% compute FOM tensors
disp(' PFOM HYDRODYNAMIC TENSORS:')
fprintf(' Assembling %d elements ...\n', nel)

fprintf('   0th order in ud:\n')
tic
T1 = myAssembly.vector_skin('Te1', 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('       1st order terms - T1: %.2f s\n',toc)

tic
Tu2 = myAssembly.matrix_skin('Teu2', 'weights', skinElements, skinElementFaces, vwater, rho);
Tudot2 = myAssembly.matrix_skin('Teudot2', 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('       2nd order terms - Tu2, Tudot2: %.2f s\n',toc)


tic;
Tuu3 = 0.5*myAssembly.tensor_skin('Teuu3',[nDOFs nDOFs nDOFs],[2 3],'weights', skinElements, skinElementFaces, vwater, rho);
Tuudot3 = myAssembly.tensor_skin('Teuudot3',[nDOFs nDOFs nDOFs],[2 3], 'weights', skinElements, skinElementFaces, vwater, rho);
Tudotudot3 = 0.5*myAssembly.tensor_skin('Teudotudot3',[nDOFs nDOFs nDOFs],[2 3], 'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('       3rd order terms - Truu3, Truudot3, Turudotudot3: %.2f s\n',toc)


fprintf('   1st order in ud:\n')

tic;
T2 = myAssembly.matrix_skin_PFOM('Te2', U,'weights', skinElements, skinElementFaces, vwater, rho);
fprintf('       2nd order terms - Tr2: %.2f s\n',toc)

tic;                            
Tu3 = myAssembly.tensor_skin_PFOM('Teu3', U, [nDOFs nDOFs nDOFs], [2 3], 'weights', skinElements, skinElementFaces, vwater, rho);
Tudot3 = myAssembly.tensor_skin_PFOM('Teudot3', U, [nDOFs nDOFs nDOFs], [2 3], 'weights', skinElements, skinElementFaces, vwater, rho);
time = toc(t0);

fprintf(' %.2f s (%.2f s)\n',toc(t0),time)
fprintf(' SPEED: %.1f el/s\n',nel/time)


% outputs
tensors.T1 = T1;
tensors.T2 = T2;
tensors.Tu2 = Tu2;
tensors.Tu3 = Tu3;
tensors.Tudot2 = Tudot2;
tensors.Tudot3 = Tudot3;
tensors.Tuu3 = Tuu3;
tensors.Tuudot3 = Tuudot3;
tensors.Tudotudot3 = Tudotudot3; 
tensors.time = time;

end


