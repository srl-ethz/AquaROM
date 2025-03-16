% cost_function
%
% Synthax:
% L =cost_function(N,tensors_hydro_PROM,u,ud,d,V)
%
% Description: Convert a constant direction vector to its reduced version
%
% INPUTS: 
% (1) N:                number of time steps             
% (2) tensors_hydro:    reduced tensors for the hydrodynamic forces
% (3) u:                solution for the nodal displacement
% (4) ud:               solution for the nodal velocities
% (5) d:                forward swimming direction vector 
% (6) V:                ROB matrix
%
%
% OUTPUTS:
% (1) Lr:   reduced cost function value    
%     
%
% Additional notes: NOT IMPLEMENTED YET !!
%
% Last modified: 19/12/2022, Mathieu Dubied, ETH Zürich

function L = cost_function(N,tensors_hydro_PROM,u,ud,d,V)
    Lr = 0;
    for i=1:N-2
        for j=1:
        eta_i = u(:,i);
        etad_i = ud(:,i);
        fhydro = double(tensors_hydro_PROM.Tr1) + ...
        double(tensors_hydro_PROM.Tru2*eta_i) + double(tensors_hydro_PROM.Trudot2*etad_i) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truu3,eta_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Truudot3,etad_i,3), eta_i,2)) + ...
        double(ttv(ttv(tensors_hydro_PROM.Trudotudot3,etad_i,3), etad_i,2));
       
        L = L - d'*fhydro;
    end
end