function h = PlotSpineAndTailElements(Nodes,Elements, ...
    SpineElementWeights,ColorSpine,TailElementWeights,ColorTail)
%--------------------------------------------------------------------------
% Purpose:
%         Plot the nominal shape (no deformation) and highlights the spine
%         elements and the tail element
% Variable Description:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           - SpineElementWeights - vector of size nElements, with 1 if the
%           element is part of the spine
%           - ColorSpine - color to use when showing the spine
%           - TailElementWeights - vector of size nElements, with 1 if the
%           element is the tail element (a single element)
%           - ColorTail - color to use when showing the tail
%
% Last modified: 10/11/2024, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------

meshcolor = 'k';
[skin,~,~,~] = getSkin3D(Elements);      
skinFaces = skin.';
nSkinFaces = size(skinFaces,1);         % total number of faces
nodePerSkinFace = size(skinFaces,2);    % number of nodes per face

X = Nodes(skinFaces',1); X = reshape(X, nodePerSkinFace, nSkinFaces);
Y = Nodes(skinFaces',2); Y = reshape(Y, nodePerSkinFace, nSkinFaces);
Z = Nodes(skinFaces',3); Z = reshape(Z, nodePerSkinFace, nSkinFaces);

view(3)
hold on;

h{1} = patch(X,Y,Z,'white','EdgeColor',meshcolor,'DisplayName','Deformed Mesh','FaceAlpha',.5);

% Spine elements __________________________________________________________
for idx=1:size(SpineElementWeights,1)
    if SpineElementWeights(idx) == 1 
        % Get the node indices for the element
        elementNodes = Elements(idx, :);
        
        % Retrieve the coordinates for the nodes of this element
        x = Nodes(elementNodes, 1);
        y = Nodes(elementNodes, 2);
        z = Nodes(elementNodes, 3);
        
        % Define the four faces of the tetrahedron
        faces = [
            1 2 3;
            1 2 4;
            1 3 4;
            2 3 4
        ];
        
        % Plot each face separately
        for j = 1:4
            patch('Vertices', [x y z], 'Faces', faces(j, :), ...
                  'FaceColor', ColorSpine, 'EdgeColor', 'black', 'FaceAlpha', 0.5);
        end
    end
end


% Tail element ____________________________________________________________
for idx=1:size(TailElementWeights,1)
    if TailElementWeights(idx) == 1 
        % Get the node indices for the element
        elementNodes = Elements(idx, :);
        
        % Retrieve the coordinates for the nodes of this element
        x = Nodes(elementNodes, 1);
        y = Nodes(elementNodes, 2);
        z = Nodes(elementNodes, 3);
        
        % Define the four faces of the tetrahedron
        faces = [
            1 2 3;
            1 2 4;
            1 3 4;
            2 3 4
        ];
        
        % Plot each face separately
        for j = 1:4
            patch('Vertices', [x y z], 'Faces', faces(j, :), ...
                  'FaceColor', ColorTail, 'EdgeColor', 'black', 'FaceAlpha', 0.5);
        end
    end
end

axis equal;
axis off;
hold off;

end
