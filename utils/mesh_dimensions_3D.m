% mesh_dimensions_3D
%
% Synthax:
% [Lx, Ly, Lz] = mesh_dimensions_3D(nodes)
%
% Description: Computes the dimensions of the FE mesh along the three axis 
%
% INPUTS: 
% (1) nodes: nodes (i.e., their position) of the FE mesh              
%
% OUTPUTS:   
% (1) Lx: dimension along the x-axis (length)
% (2) Ly: dimension along the y-axis (depth)
% (3) Lz: dimension along the z-axis (height)
%
% Additional notes: -
%
% Last modified: 02/09/2024, Mathieu Dubied, ETH Zürich
function [Lx, Ly, Lz] = mesh_dimensions_3D(nodes)
    Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % dimension along the x-axis (length)
    Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % dimension along the y-axis (depth)
    Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % dimension along the z-axis (height)
end