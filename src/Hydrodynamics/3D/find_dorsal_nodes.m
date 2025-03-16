% find_dorsal_nodes
%
% Synthax:
% [allDorsalNodesIdx,matchedDorsalNodesIdx,dorsalNodesElementsVec,matchedDorsalNodesZPos] = ...
%   find_dorsal_nodes(elements, nodes, spineElements, nodeIdxPosInElements)
%
% Description: 
% Find the dorsal nodes matching the spine nodes/elements
%
% INPUTS
%   - elements:             a 2D array in which each row contains the
%                           nodes' indexes of an element
%   - nodes:                a 2D array containing the x,y and z positions 
%                           of each node
%   - spineElements:        array containing the indexes of the spine
%                           elements
%   - nodeIdxPosInElements: matrix of size nElements x 2. For the row
%                           corresponding to spine elements, it gives the
%                           column index (1 to 4) in the element of the 
%                           node close to the tail (1st column) and the 
%                           node close to the head (2nd column).
%
% OUTPUT:
%   - allDorsalNodesIdx:    array containing the indexes of the dorsal
%                           nodes (all of them)
%   - matchedDorsalNodesIdx:array containing the indexes of the dorsal
%                           nodes that have been matched to a spine
%                           element
%   - dorsalNodesElementsVec: array of size nElements, in which the entry
%                           corresponding to spine elements is filled with
%                           the matched dorsal node index.
%   - matchedDorsalNodesZPos: array of size nElements, in which the entry
%                           corresponding to spine elements is filled with
%                           the matched dorsal node position. For each
%                           spine element, the matching with a (single)
%                           dorsal node is performed considering the spine 
%                           nodes located toward the tail end of the fish.

%
% Additional notes: Only test for candidate dorsal nodes that are aligned
% with the spine nodes. 
%
% Last modified: 02/05/2024, Mathieu Dubied, ETH Zurich
function [allDorsalNodesIdx,matchedDorsalNodesIdx,dorsalNodesElementsVec,matchedDorsalNodesZPos] = ....
    find_dorsal_nodes(elements, nodes, spineElements, nodeIdxPosInElements)
       
    % GET ALL DORSAL NODES ________________________________________________
    
    % select nodes located in the upper half xz plane (y=0)
    upHalfIdx = find(nodes(:,2)==0 & nodes(:,3)>0);

    % get nodes that are part of the skin
    [skin,~,~,~] = getSkin3D(elements);
    
    % find common node indexes between skin and up Half Idx
    % note: it also contains the skin nodes on the tail and head
    allDorsalNodesIdx = intersect(upHalfIdx,skin);
  
    % MATCHING SPINE ELMENTS WITH A SINGLE DORSAL NODES ___________________
    nElements = length(elements(:,1));
    nSpineEl = length(spineElements);
    matchedDorsalNodesIdx = zeros(nSpineEl,1);
    matchedDorsalNodesZPos = zeros(nElements,1);
    dorsalNodesElementsVec = zeros(nElements,1);
    for i=1:nSpineEl
        % get spine node position
        el = spineElements(i);
        spineNodeIdx = elements(el,nodeIdxPosInElements(el,1));
        spineNodePos = nodes(spineNodeIdx,1);   % x position
        
        % get x and z position of dorsal nodes
        nodes2searchX = nodes(allDorsalNodesIdx,1);
        nodes2searchZ = nodes(allDorsalNodesIdx,3);

        % match spine node with dorsal node (min x distance)
        xTargetVec = repmat(spineNodePos,size(nodes2searchX,1),1);
        sn = abs(nodes2searchX - xTargetVec) ;
        candidatesIdx = find(sn(:)==min(sn));
        if length(candidatesIdx) == 1
            ind = candidatesIdx;
        else
            % multiple candidates for a single spine node, happen at the
            % tail of the fish. Select the one with the largest z value
            zPosOfCandidates = nodes2searchZ(candidatesIdx);
            [~,idxInZPosList] = max(zPosOfCandidates);
            ind = candidatesIdx(idxInZPosList);
        end

        % get matched dorsal node index
        matchedDorsalNodesIdx(i) = allDorsalNodesIdx(ind);
        matchedDorsalNodesZPos(el) = nodes(allDorsalNodesIdx(ind),3);
        dorsalNodesElementsVec(el) = allDorsalNodesIdx(ind);
          
    end
  
end