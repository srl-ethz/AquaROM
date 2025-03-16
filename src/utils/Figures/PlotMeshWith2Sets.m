function h = PlotMeshWith2Sets(Nodes,Elements, ...
    ElementSet1,ColorSet1,ElementSet2,ColorSet2, varargin)
%--------------------------------------------------------------------------
% Purpose:
%           To plot 2D and 3D Finite Element Method Mesh and highlight 2
%           sets of elements in two different colors
%         
% Variable Description:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           - ElementSet1 - vector of size nElements, with 1 if the
%           element is part of the set 1
%           - ColorSet1 - color to use when showing the first element set
%           - ElementSet2 - vector of size nElements, with 1 if the
%           element is part of the set 2
%           - ColorSet2 - color to use when showing the second element set
%
% Additional notes: only tested for QUAD4, TRI3, and TET4 elements. Based
% on the code PlotMesh. TO DO test for TET4
%
% Last modified: 17/12/2024, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------

[alphaSet1, alphaSet2, show] = parse_inputs(varargin{:});

dimension = size(Nodes,2) ;  % Dimension of the mesh
nel = size(Elements,1) ;                  % number of elements
nnode = length(Nodes) ;          % total number of nodes in system
nnel = size(Elements,2);                % number of nodes per element

% 3D PLOTS ________________________________________________________________
if dimension == 3   % For 3D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:3)));
    
    
    if elementdim == 3 % solid in 3D when we simply plot the skin elements
        faces = getSkin3D(Elements);
        Elements = faces.';
        nel = size(Elements,1);      % total number of faces
        nnel = size(Elements,2);     % number of nodes per face
    end
    
    X = Nodes(Elements',1); X = reshape(X, nnel, nel);
    Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
    Z = Nodes(Elements',3); Z = reshape(Z, nnel, nel);
    
    if elementdim ~= 3 % surface or line in 3D        
        [X,Y,Z] = tune_coordinates(X,Y,Z);
    end
    
    p = patch(X,Y,Z,'w','FaceAlpha',1.0,'EdgeAlpha',1,...
        'EdgeColor','k','LineStyle','-','DisplayName','Mesh');
    view(3)
    
% 2D PLOTS ________________________________________________________________
elseif dimension == 2           % For 2D plots
    elementdim = rank(diff(Nodes(Elements(1,:),1:2))); % dimension of element
    hold on
    if elementdim == 2
        X = Nodes(Elements',1); X = reshape(X, nnel, nel);
        Y = Nodes(Elements',2); Y = reshape(Y, nnel, nel);
        p = patch(X,Y,'w','DisplayName','Mesh','LineWidth',0.01,'EdgeAlpha',0.2); 
    else % line
        p = plot(Nodes(:,1),Nodes(:,2),'.-k', 'Markersize',10);% no linewidth
    end
    
    % element set 1
    for idx=1:size(ElementSet1,1)
        if ElementSet1(idx) == 1 
            nodeIdx = Elements(idx,:);
            p = patch(Nodes(nodeIdx,1),Nodes(nodeIdx,2),ColorSet1, 'FaceAlpha',alphaSet1,'LineWidth',0.01,'EdgeAlpha',0.2);
        end
    end
    
    % element set 2
    for idx=1:size(ElementSet2,1)
        if ElementSet2(idx) == 1 
            nodeIdx = Elements(idx,:);
            p = patch(Nodes(nodeIdx,1),Nodes(nodeIdx,2),ColorSet2, 'FaceAlpha',alphaSet2,'LineWidth',0.01,'EdgeAlpha',0.2);
        end
    end
    
    
end
% display Node numbers and Element numbers
if show ~= 0
    k = 1:nnode ;
    nd = k' ;
    if size(Nodes,2)==2
        for i = 1:nel
            text(X(:,i),Y(:,i),int2str(nd(Elements(i,:))),'fontsize',8,'color','k');
            text(mean(X(:,i)),mean(Y(:,i)),int2str(i),'fontsize',10,'color','r') ;
        end
    elseif size(Nodes,2)==3
        for i = 1:nel
            text(X(:,i),Y(:,i),Z(:,i),int2str(nd(Elements(i,:))),'fontsize',8,'color','k');
            text(mean(X(:,i)),mean(Y(:,i)),mean(Z(:,i)),int2str(i),'fontsize',10,'color','r') ;
        end
        p.FaceAlpha = 0.2;
    end
end
% set(gca,'XTick',[]) ; set(gca,'YTick',[]); set(gca,'ZTick',[]);
rotate3d on;
axis equal;
axis off;
end

function [X,Y,Z] = tune_coordinates(X,Y,Z)
% ONLY for UNDEFORMED mesh:
% to avoid visualization problems when the surface is contained in
% a plane parallel to the main planes (xy,xz,yz):
if sum(X(:)==X(1))==numel(X)
    X(1)=X(1)+1e-15;
elseif sum(Y(:)==Y(1))==numel(Y)
    Y(1)=Y(1)+1e-15;
elseif sum(Z(:)==Z(1))==numel(Z)
    Z(1)=Z(1)+1e-15;
end
end

%% PARSING INPUTS (VARARGIN) ______________________________________________
function [alphaSet1, alphaSet2, show] = parse_inputs(varargin)

p = inputParser;
addParameter(p,'alphaSet1', 0.6);
addParameter(p,'alphaSet2', 0.6);
addParameter(p,'show', 0);

parse(p,varargin{:});

alphaSet1 = p.Results.alphaSet1;
alphaSet2= p.Results.alphaSet2;
show = p.Results.show;

end