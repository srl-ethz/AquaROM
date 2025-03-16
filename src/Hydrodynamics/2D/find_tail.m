% find_tail
%
% Synthax:
% [tailNodeOfTailElement, tailElementIdx, tailElementWeights] = find_tail_(elements, nodes, spineElements,nodeIdxPosInElements)
%
% Description: 
% Find the tail node and element, return related weights vector
%
% INPUTS
%   - elements: a 2D array in which each row contains the nodes' indexes of
%               an element
%   - nodes:    a 2D array containing the x and y positions of each node
%
%   - spineElements         array containing the indexes of the spine elements.
%                           In case there are two elements are sharing a face 
%                           on the spine, only the upper element is considered.
%                           This allows to avoid applying the force twice.
%   - nodeIdxPosInElements: matrix of size nElements x 2. For the row
%                           corresponding to spine elements, it gives the
%                           column index (1 to 3/4) in the element of the 
%                           node close to the tail (1st column) and the 
%                           node close to the head (2nd column).
% OUTPUTS
%   - tailNodeOfTailElement:tail node of the tail element, so the tail the
%                           node at the tail of the fish
%   - tailElementIdx:       index of the element being in the spine and at
%                           the tail of the fish
%   - tailElementWeights:   array of length nElements, with 1 if element is
%                           the tail element, 0 else 
%
% Additional notes: -
%
% Last modified: 06/10/2023, Mathieu Dubied, ETH Zurich

function [tailNodeOfTailElement, tailElementIdx, tailElementWeights] = find_tail(elements, nodes, spineElements,nodeIdxPosInElements)
    nElements = length(elements(:,1));
    tailElementWeights = zeros(nElements,1);

    % get useful variables
    nodesOfSpineElements = elements(spineElements,:);
    tailNodesInSpineElements = nodeIdxPosInElements(spineElements,1);
    % linear indexing to get tail nodes
    I = (1:size(nodesOfSpineElements, 1)).';
    k = sub2ind(size(nodesOfSpineElements), I, tailNodesInSpineElements);
    tailNodes = nodesOfSpineElements(k);
    
    % find nodes with lowest x position among tail nodes
    [~,xMinIndex] = min(nodes(tailNodes,1));
    
    % find tail element
    tailElementIdx = spineElements(xMinIndex);
    tailNodeOfTailElement = elements(tailElementIdx,tailNodesInSpineElements(xMinIndex));
    
    % weigth output
    tailElementWeights(tailElementIdx) = 1;

end
