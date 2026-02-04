function h = HighlightElement(Nodes,Elements,ElementWeights,Color)
%--------------------------------------------------------------------------
% Purpose:
%         Plot the nominal shape (no deformation) and highlights the spine
%         elements and the tail element
% Variable Description:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           - ElementWeights - vector of size nElements, with 1 if the
%           element is to be highlighted
%           - Color - color to use when highlighting the element(s)
%
% Last modified: 03/02/2026, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------

meshcolor = 'k';
[skin,~] = getSkin3D(Elements);      
skinFaces = skin.';
nSkinFaces = size(skinFaces,1);         % total number of faces
nodePerSkinFace = size(skinFaces,2);    % number of nodes per face

X = Nodes(skinFaces',1); X = reshape(X, nodePerSkinFace, nSkinFaces);
Y = Nodes(skinFaces',2); Y = reshape(Y, nodePerSkinFace, nSkinFaces);
Z = Nodes(skinFaces',3); Z = reshape(Z, nodePerSkinFace, nSkinFaces);

view(3)
hold on;

h{1} = patch(X,Y,Z,'white','EdgeColor',meshcolor,'DisplayName','Deformed Mesh','FaceAlpha',.5);

% Element to highlight ____________________________________________________
for idx=1:size(ElementWeights,1)
    if ElementWeights(idx) == 1 
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
                  'FaceColor', Color, 'EdgeColor', 'black', 'FaceAlpha', 0.5);
        end
    end
end

axis equal;
axis off;
hold off;

end
