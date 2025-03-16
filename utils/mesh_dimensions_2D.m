% mesh_dimensions_2D
%
% Synthax:
% [Lx, Ly] = mesh_dimensions_2D(nodes)
%
% Description: Computes the dimensions of the FE mesh along the two axis 
%
% INPUTS: 
% (1) nodes: nodes (i.e., their position) of the FE mesh              
%
% OUTPUTS:   
% (1) Lx: dimension along the x-axis (length)
% (2) Ly: dimension along the y-axis (height)
%
% Additional notes: -
%
% Last modified: 17/012/2024, Mathieu Dubied, ETH Zürich
function [Lx, Ly] = mesh_dimensions_2D(nodes)
    Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % dimension along the x-axis (length)
    Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % dimension along the y-axis (height)
end