% shape_variations_3D
%
% Synthax:
% [thinFish,shortFish,linearTail,longTail,shortTail,linearHead,longHead,shortHead] 
% = shape_variations_3D(nodes,Lx,Ly)
%
% Description: Defines the shape variations used for the 3D case study
%
% INPUTS:              
% (1) nodes:    nodes and their coordinates
% (2) Lx:       size of the nominal shape in the x-direction
% (3) Ly:       size of the nominal shape in the y-direction
% (4) Lz:       size of the nominal shape in the z-direction
%
% OUTPUTS:
% (1) [...]:    11 shape variations
%     
%
% Last modified: 06/02/2025, Mathieu Dubied, ETH Zurich
function [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch, ...
    y_tail,y_head,y_linLongTail,y_ellipseFish,xz_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz)

    % (1) thinner fish (y direction)
    nodes_projected = [nodes(:,1), nodes(:,2)*0, nodes(:,3)];   % projection on y-axis
    yDif = nodes_projected(:,2) - nodes(:,2);                   % y-difference projection vs nominal
    y_thinFish = zeros(numel(nodes),1);                         
    y_thinFish(2:3:end) = yDif;                                 
    
    % (2) smaller fish (z direction)
    nodes_projected = [nodes(:,1), nodes(:,2), nodes(:,3)*0];   % projection on z-axis
    zDif = nodes_projected(:,3) - nodes(:,3);                   % z-difference projection vs nominal
    z_smallFish = zeros(numel(nodes),1);            
    z_smallFish(3:3:end) = zDif;                   
    
    % (3) smaller fish tail (z direction)
    elCenter = [-Lx*0.7;0];
    a = Lx*0.3;
    b = Lz*0.5;
    nodes_ellipse = nodes;
    zDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1) 
        if nodes(n,1) <= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), nodes(n,2), sign(nodes(n,3)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
            zDif(n) = (nodes_ellipse(n,3) - max(nodes(:,3)).*sign(nodes(n,3))).*abs(nodes(n,3))/(Lz*0.5);                   % z-difference projection vs nominal
        end
    end
    z_tail = zeros(numel(nodes),1);
    z_tail(3:3:end) = real(zDif);
    
    
    % (4) smaller fish head (z direction)
    elCenter = [-Lx*0.3;0];
    a = Lx*0.3;
    b = Lz*0.5;
    nodes_ellipse = nodes;
    zDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), nodes(n,2), sign(nodes(n,3)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2)];    % projection on ellipse
            zDif(n) = (nodes_ellipse(n,3) - max(nodes(:,3)).*sign(nodes(n,3))).*abs(nodes(n,3))/(Lz*0.5);                   % z-difference projection vs nominal
        end
    end
    z_head = zeros(numel(nodes),1);
    z_head(3:3:end) = real(zDif);
    
    % (5) smaller linear long fish tail (z direction)
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
    z_linLongTail = zeros(numel(nodes),1);
    z_linLongTail(3:3:end) = real(zDif);
    
    % (6) notch before the tail (z direction)
    xFront = -0.3*Lx;
    xBack = -0.9*Lx;
    xNotch = -Lx*0.7;
    slopeFront = 0.5*Lz/(xFront-xNotch);
    slopeBack = -0.5*Lz/(xBack-xNotch);
    interceptFront = 0.5*Lz-xFront*slopeFront;
    interceptBack = 0.5*Lz-xBack*slopeBack;
    nodes_proj = nodes;
    zDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1) 
        if nodes(n,1) <= xNotch && nodes(n,1)> xBack
            nodes_proj(n,:) = [nodes(n,1), nodes(n,2), slopeBack*nodes(n,1)+interceptBack];                         % projection on linear curve
            zDif(n) = abs(nodes_proj(n,3) - max(nodes(:,3))).*-sign(nodes(n,3)).*abs(nodes(n,3))/(Lz*0.5);  % z-difference projection vs nominal
        elseif nodes(n,1) > xNotch && nodes(n,1)< xFront
            nodes_proj(n,:) = [nodes(n,1), nodes(n,2), slopeFront*nodes(n,1)+interceptFront];                         % projection on linear curve
            zDif(n) = abs(nodes_proj(n,3) - max(nodes(:,3))).*-sign(nodes(n,3)).*abs(nodes(n,3))/(Lz*0.5);  % z-difference projection vs nominal
        end 
    end
    z_notch = zeros(numel(nodes),1);
    z_notch(3:3:end) = real(zDif);

    % (7) smaller fish tail (y direction)
    elCenter = [-Lx*0.7;0];
    a = Lx*0.3;
    b = Ly*0.5;
    nodes_ellipse = nodes;
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2), nodes(n,3)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);                   % y-difference projection vs nominal
        end
    end
    y_tail = zeros(numel(nodes),1);
    y_tail(2:3:end) = real(yDif);
    
    % (8) smaller fish head (y direction)
    elCenter = [-Lx*0.3;0];
    a = Lx*0.3;
    b = Ly*0.5;
    nodes_ellipse = nodes;
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) >= elCenter(1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,1)-elCenter(1)).^2), nodes(n,3)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);                   % y-difference projection vs nominal
        end
    end
    y_head = zeros(numel(nodes),1);
    y_head(2:3:end) = real(yDif);
       
    % (9) smaller linear long fish tail(y direction)
    xStart = -Lx*0.3;
    slope = 0.5*Ly/(0.7*Lx);
    intercept = 0.5*Ly+0.3*Lx*slope;
    nodes_proj = nodes;
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
        if nodes(n,1) <= xStart
            nodes_proj(n,:) = [nodes(n,1), slope*nodes(n,1)+intercept, nodes(n,3)];                             % projection on linear curve
            yDif(n) = abs(nodes_proj(n,2) - max(nodes(:,2))).*-sign(nodes(n,2)).*abs(nodes(n,2))/(Ly*0.5);      % y-difference projection vs nominal
        end
    end
    y_linLongTail = zeros(numel(nodes),1);
    y_linLongTail(2:3:end) = real(yDif);
    
    % (10) ellipse in yz plane over whole fish
    elCenter = [0;0];
    a = Lz*0.5;
    b = Ly*0.5;
    nodes_ellipse = nodes;
    yDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1)
            nodes_ellipse(n,:) = [nodes(n,1), sign(nodes(n,2)).*b/a.*sqrt(a^2-(nodes(n,3)-elCenter(1)).^2), nodes(n,3)];    % projection on ellipse
            yDif(n) = (nodes_ellipse(n,2) - max(nodes(:,2)).*sign(nodes(n,2))).*abs(nodes(n,2))/(Ly*0.5);                   % y-difference projection vs nominal
    end
    y_ellipseFish = zeros(numel(nodes),1);
    y_ellipseFish(2:3:end) = real(yDif);
    
    % (11) concave tail
    xDif = zeros(size(nodes,1),1);
    zDif = zeros(size(nodes,1),1);
    for n = 1:size(nodes,1) 
        if nodes(n,1) == -Lx
            xDif(n) = -1*abs(nodes(n,3));  
            zDif(n) = 0.4*nodes(n,3); 
        end
    end
    xz_concaveTail = zeros(numel(nodes),1);
    xz_concaveTail(1:3:end) = xDif;
    xz_concaveTail(3:3:end) = zDif;
    
end








