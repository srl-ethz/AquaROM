% hydro_force_TRI3
%
% Synthax:
% force = hydro_force_TRI3(myAssembly, elements, skinElements, skinElementFaces, vwater, rho)
%
% Description: 
% This function computes the hydrodynamic forces at the Assembly level,
% using the original nonlinear formulation (not in form of a polynomial).
%
%
% INPUTS
%   - myAssembly: Assembly from YetAnotherFEcode.
%   - skinElements: logical vector (0 or 1), with 1 at position i
%                   indicating that the i-th element is part of the skin
%   - skinElementFaces: each row of this matrix show which face of the
%                   element is part of the skin (1: nodes 1 and 2, 2: nodes
%                   2 and 3, 3: nodes 3 and 1)
%   - vwater: 2D vector for the water velocity
%   - rho: rho of the fluid (water)
%
% OUTPUT:
%   force: a struct variable with the following fields:
%       .f              hydrodynamic force at the Assembly level        
%      	.time           computational time
%     
%
% Additional notes:
%   - List of currently supported elements: 
%     TRI3
%
% Last modified: 16/04/2023, Mathieu Dubied, ETH Zurich
function force = hydro_force_TRI3(Assembly, skinElements, skinElementFaces, vwater, rho, c, q, qd)
    u = Assembly.unconstrain_vector(q);
    ud = Assembly.unconstrain_vector(qd);
    force = Assembly.vector_hydro('drag_force_full', 'weights', skinElements, skinElementFaces, vwater, rho,c,u,ud) + ...
            Assembly.vector_hydro('thrust_force_full', 'weights', skinElements, skinElementFaces, vwater, rho,c,u,ud);
    force = Assembly.constrain_vector(force);
end

