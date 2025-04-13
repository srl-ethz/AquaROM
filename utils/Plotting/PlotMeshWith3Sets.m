function h = PlotMeshWith3Sets(Nodes,Elements, ...
    ElementSet1,ColorSet1,ElementSet2,ColorSet2,ElementSet3,ColorSet3)
%--------------------------------------------------------------------------
% Purpose:
%         Plot the mesh and highlights the 3 different sets of elements
% Variable Description:
%           Nodes - The nodal coordinates of the mesh
%           -----> Nodes = [X Y Z]
%           Elements - The nodal connectivity of the elements
%           -----> Elements = [node1 node2......]
%           - ElementSet1 - vector of size nElements, with 1 if the
%           element is part of the first muscle
%           - ColorSet1 - color to use when showing the first muscle
%           - ElementSet2 - vector of size nElements, with 1 if the
%           element is part of the second muscle
%           - ColorSet2 - color to use when showing the second muscle
%           - ElementSet3 - vector of size nElements, with 1 if the
%           element is part of the elements that are constrained
%           - ColorSet3 - color to use when showing the constrained
%           elements
%
% Additional notes: only tested for TET4 elements
%
% Last modified: 16/03/2025, Mathieu Dubied, ETH Zurich
%--------------------------------------------------------------------------

meshcolor = 'k';
alphaValue = 0.7;

originalElements = Elements;
[skin,~,skinElements,skinElementFaces] = getSkin3D_ext(Elements);      
skinFaces = skin.';
nSkinFaces = size(skinFaces,1);         % total number of faces
nodePerSkinFace = size(skinFaces,2);    % number of nodes per face

X = Nodes(skinFaces',1); X = reshape(X, nodePerSkinFace, nSkinFaces);
Y = Nodes(skinFaces',2); Y = reshape(Y, nodePerSkinFace, nSkinFaces);
Z = Nodes(skinFaces',3); Z = reshape(Z, nodePerSkinFace, nSkinFaces);

% open View->Camera Toolbar, and View->Properties Inspector (Axes,
% view(-1.5078,39.6338)   
% camproj('orthographic')
% campos([-0.949358102339665,-1.789111161270594,1.491057614878067])
% camtarget([-0.160061975241042,-0.004309698762181,0.011789814756804])
% camup([0.244316050410686,0.552461400881531,0.79692914870002])
% camva(11.0954365)
view(3)
hold on;

h{1} = patch(X,Y,Z,'white','EdgeColor',meshcolor,'DisplayName','Deformed Mesh');
hIdx = 2;

% First muscle ____________________________________________________________
for idx=1:size(ElementSet1,1)
    % check if an element is part of the skin elements and actuation
    % elements

    if ElementSet1(idx) == 1 && skinElements(idx) == 1
        % find nodes of the skin faces of this element
        faceNodes = faceFromsSkinElement(originalElements(idx,:),skinElementFaces(idx,:));
        % find idx in `skin' (i.e. skin faces) that corresponds to the
        % faceNodes we just found (order can be different)
        [i1, ~]=find(skinFaces==faceNodes(1,1));
        [i2, ~]=find(skinFaces==faceNodes(1,2));
        [i3, ~]=find(skinFaces==faceNodes(1,3));
        i12 = intersect(i1,i2);
        ColumnIdx = intersect(i12,i3);

        h{hIdx} = patch(X(ColumnIdx*3-2:ColumnIdx*3),Y(ColumnIdx*3-2:ColumnIdx*3),Z(ColumnIdx*3-2:ColumnIdx*3),ColorSet1,'EdgeColor',meshcolor);
        set(h{hIdx},'FaceAlpha',alphaValue)
        hIdx = hIdx + 1;

        % possibly a second element in this face being par of the skin
        if faceNodes(2,1) ~= 0
            [i1, ~]=find(skinFaces==faceNodes(2,1));
            [i2, ~]=find(skinFaces==faceNodes(2,2));
            [i3, ~]=find(skinFaces==faceNodes(2,3));
            i12 = intersect(i1,i2);
            ColumnIdx = intersect(i12,i3);

            h{hIdx} = patch(X(ColumnIdx*3-2:ColumnIdx*3),Y(ColumnIdx*3-2:ColumnIdx*3),Z(ColumnIdx*3-2:ColumnIdx*3),ColorSet1,'EdgeColor',meshcolor);
            set(h{hIdx},'FaceAlpha',alphaValue)
            hIdx = hIdx + 1;
        end

    end
end

% Second muscle ___________________________________________________________
for idx=1:size(ElementSet2,1)
    % check if an element is part of the skin elements and actuation
    % elements

    if ElementSet2(idx) == 1 && skinElements(idx) == 1
        % find nodes of the skin faces of this element
        faceNodes = faceFromsSkinElement(originalElements(idx,:),skinElementFaces(idx,:));
        % find idx in `skin' (i.e. skin faces) that corresponds to the
        % faceNodes we just found (order can be different)
        [i1, ~]=find(skinFaces==faceNodes(1,1));
        [i2, ~]=find(skinFaces==faceNodes(1,2));
        [i3, ~]=find(skinFaces==faceNodes(1,3));
        i12 = intersect(i1,i2);
        ColumnIdx = intersect(i12,i3);

        h{hIdx} = patch(X(ColumnIdx*3-2:ColumnIdx*3),Y(ColumnIdx*3-2:ColumnIdx*3),Z(ColumnIdx*3-2:ColumnIdx*3),ColorSet2,'EdgeColor',meshcolor);
        set(h{hIdx},'FaceAlpha',alphaValue)
        hIdx = hIdx + 1;

        % possibly a second element in this face being par of the skin
        if faceNodes(2,1) ~= 0
            [i1, ~]=find(skinFaces==faceNodes(2,1));
            [i2, ~]=find(skinFaces==faceNodes(2,2));
            [i3, ~]=find(skinFaces==faceNodes(2,3));
            i12 = intersect(i1,i2);
            ColumnIdx = intersect(i12,i3);

            h{hIdx} = patch(X(ColumnIdx*3-2:ColumnIdx*3),Y(ColumnIdx*3-2:ColumnIdx*3),Z(ColumnIdx*3-2:ColumnIdx*3),ColorSet2,'EdgeColor',meshcolor);
            set(h{hIdx},'FaceAlpha',alphaValue)
            hIdx = hIdx + 1;
        end

    end
end

% Constrained elements ____________________________________________________
for idx=1:size(ElementSet3,1)
    % check if an element is part of the skin elements and actuation
    % elements

    if ElementSet3(idx) == 1 && skinElements(idx) == 1
        % find nodes of the skin faces of this element
        faceNodes = faceFromsSkinElement(originalElements(idx,:),skinElementFaces(idx,:));
        % find idx in `skin' (i.e. skin faces) that corresponds to the
        % faceNodes we just found (order can be different)
        [i1, ~]=find(skinFaces==faceNodes(1,1));
        [i2, ~]=find(skinFaces==faceNodes(1,2));
        [i3, ~]=find(skinFaces==faceNodes(1,3));
        i12 = intersect(i1,i2);
        ColumnIdx = intersect(i12,i3);

        h{hIdx} = patch(X(ColumnIdx*3-2:ColumnIdx*3),Y(ColumnIdx*3-2:ColumnIdx*3),Z(ColumnIdx*3-2:ColumnIdx*3),ColorSet3,'EdgeColor',meshcolor);
        set(h{hIdx},'FaceAlpha',alphaValue)
        hIdx = hIdx + 1;

        % possibly a second element in this face being par of the skin
        if faceNodes(2,1) ~= 0
            [i1, ~]=find(skinFaces==faceNodes(2,1));
            [i2, ~]=find(skinFaces==faceNodes(2,2));
            [i3, ~]=find(skinFaces==faceNodes(2,3));
            i12 = intersect(i1,i2);
            ColumnIdx = intersect(i12,i3);

            h{hIdx} = patch(X(ColumnIdx*3-2:ColumnIdx*3),Y(ColumnIdx*3-2:ColumnIdx*3),Z(ColumnIdx*3-2:ColumnIdx*3),ColorSet3,'EdgeColor',meshcolor);
            set(h{hIdx},'FaceAlpha',alphaValue)
            hIdx = hIdx + 1;

        end

    end
end

axis equal;
axis off;
hold off;

end

function nodes = faceFromsSkinElement(element,faceNumber)
    nodes = [0 0 0; 0 0 0];
    node = faceNumber(1);
    nodes(1,1) = element(node);
    for i = 1:2
        node = node + 1;
        if node == 5
            node = 1;
        end
        nodes(1,i+1) = element(node);
    end
    
    if faceNumber(2) ~= 0
        node = faceNumber(2);
        nodes(2,1) = element(node);
        for i = 1:2
            node = node + 1;
            if node == 5
                node = 1;
            end
            nodes(2,i+1) = element(node);
        end
    end
end
