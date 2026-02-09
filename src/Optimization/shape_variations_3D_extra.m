% shape_variations_3D_extra
%
% Synthax:
% [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch, ...
%    y_tail,y_head,y_linLongTail,y_ellipseFish,xz_concaveTail]
%       = shape_variations_3D_extra(nodes,Lx,Ly)
%
% Description: Defines a few extra shape variations
%
% INPUTS:              
% (1) nodes:    nodes and their coordinates
% (2) Lx:       size of the nominal shape in the x-direction
% (3) Ly:       size of the nominal shape in the y-direction
% (4) Lz:       size of the nominal shape in the z-direction
%
% OUTPUTS:
% (1) [...]:    13 shape variations
%     
%
% Last modified: 09/02/2026, Mathieu Dubied, ETH Zurich
function [z_tail_1, z_tail_2] = ...
    shape_variations_3D_extra(nodes,Lx,Ly,Lz)
                               
    
    % (1) z_tail_1 (z direction, bulky tail)
    x_limit = - 0.5*Lx;
    nodes_projected = nodes;
    for n = 1:size(nodes,1)
        if nodes(n,1) <= x_limit
            nodes_projected(n,:) = [nodes(n,1), nodes(n,2), nodes(n,3)*0];  % projection on z-axis
        end
    end
    zDif = nodes_projected(:,3) - nodes(:,3);   % z-difference projection vs nominal
    z_tail_1 = zeros(numel(nodes),1);            
    z_tail_1(3:3:end) = zDif;                   
    
    % (2) z_tail_2 (z direction, linear tail)
    xStart = -Lx*0.3;
    slope = 0.5*Lz/(0.7*Lx);
    intercept = 0.5*Lz+0.3*Lx*slope;
    nodes_proj = nodes;
    zDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1) 
        if nodes(n,1) <= xStart
            nodes_proj(n,:) = [nodes(n,1), nodes(n,2), slope*nodes(n,1)+intercept];                         % projection on linear curve
            zDif(n) = abs(nodes_proj(n,3) - max(nodes(:,3))).*-sign(nodes(n,3)).*abs(nodes(n,3))/(Lz*0.5);  % z-difference projection vs nominal
        end
    end
    z_tail_2 = zeros(numel(nodes),1);
    z_tail_2(3:3:end) = real(zDif);
   
    
end








