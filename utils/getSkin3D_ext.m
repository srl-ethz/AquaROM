% -------------------------------------------------------------------------
% getSkin3D(elements): computes the external faces of the mesh. 
% INPUTS:
%   elements: table of elements (Nelements X Ndofs)
% OUTPUTS:
%   skin: table of nodes labels of external faces. Each column contains the
%         nodes of one face.
%   allfaces: as skin, but with all the faces.
%   skinElements: logical vector of size (nElements X 1), 1 meaning the element
%                 has a face which is part of the skin.  
%   skinElementFaces: describes which face(s) of the element is part of the
%                     skin. The faces of a given element are numerated the
%                     same way they are ordered. A single element can have
%                     up to 2 skin faces.
%
% Supported elements: TET4, maybe more as adapted from the master branch
% version
% Note: a large portion of the code is taken from getSkin3D from the master
% branch
% Last modified: 24/08/2022, Mathieu Dubied, ETH Zurich
% -------------------------------------------------------------------------
function [skin,allfaces,skinElements,skinElementFaces] = getSkin3D_ext(elements)

nnel  = size(elements,2); % number of nodes per element

% select faces order
switch nnel
    case {4,10} % tetrahedra
        faces = [1 2 4; 2 3 4; 1 3 4; 1 2 3];
    case 15 % wedges
        faces = [4 10 5 11 6 12; 1 7 2 8 3 9; 1 7 2 5 10 4; 2 8 3 6 11 5; 1 9 3 6 12 4];
    case {8,20} % hexahedra
        faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
end

% build face matrix
N = size(elements,1);
FACES = zeros(size(faces,2),N*size(faces,1));
count = 1;
for ii = 1:N
    vertici = elements(ii,:);
    for jj = 1:size(faces,1)
        FACES(:,count) = vertici(faces(jj,:))';
        count = count+1;
    end
end

% function [x1,icc] = OnlyNotRepeatedColumns(x)____________________________
x = FACES;
xs = sort(x,1);
[~,ia,ic] = unique(xs','rows','stable');
ia = ia';
ic = ic';
ic = ia(ic); % a unique index is given to each unique column
icc = hist(ic,1:length(ic));	% counts the occurences 
                                % (occ. after the first appear as zeros)
icc = icc == 1;                 % take only terms occurring ONCE
% x1 = xs(:,icc);__________________________________________________________

% remove all the repeated faces (internal), so to plot only external ones
skin = FACES(:,icc);
allfaces = FACES;


% find elements with faces being part of the skin
skinElements = zeros(N,1);
skinElementFaces = zeros(N,3);
indexSkinElements = 1;
skinMembers = ismember(allfaces.',skin','rows');

for ii = 1:size(skinMembers,1)
    if skinMembers(ii) == 1
        skinElements(indexSkinElements) = 1;
        fN = faceNumber(elements(indexSkinElements,:), allfaces(:,ii));
        if skinElementFaces(indexSkinElements,1) == 0
            skinElementFaces(indexSkinElements,1) = fN;
        else
            if skinElementFaces(indexSkinElements,2) == 0
                skinElementFaces(indexSkinElements,2) = fN;
            else
                skinElementFaces(indexSkinElements,3) = fN;
            end
            
        end
    end 

    if mod(ii,size(faces,1)) == 0
        indexSkinElements = indexSkinElements + 1;
    end
end 

% nested function faceNumber
% INPUTS: -element, vector containing the nodes' numbers of `element'
%         -nodes, vector containing the nodes in which we are interested
% OUTPUT: the face's number of `element' that corresponds to `nodes'
    function fN = faceNumber(element,nodes)
        fN = 0;
        for i=1:size(element,2)
            if i == size(element,2)
                next=1;
            else
                next = i+1;
            end
            if i+1 == size(element,2)
                nextnext = 1;
            elseif i == size(element,2)
                nextnext = 2;
            else
                nextnext = i+2;
            end
     

            if isequal([element(i);element(next);element(nextnext)],nodes) || ...
               isequal([element(i);element(nextnext);element(next)],nodes) || ...
               isequal([element(next);element(i);element(nextnext)],nodes) || ...
               isequal([element(next);element(nextnext);element(i)],nodes) || ...
               isequal([element(nextnext);element(i);element(next)],nodes) || ...
               isequal([element(nextnext);element(next);element(i)],nodes)
                fN = i;
            end           
        end
    end



end 

