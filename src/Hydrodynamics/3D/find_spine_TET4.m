% find_spine_TET4
%
% Synthax:
% [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TET4(elements, nodes)
%
% Description: 
% Find the spine nodes and elements, return related weights vector
%
% INPUTS
%   - elements: a 2D array in which each row contains the nodes' indexes of
%               an element
%   - nodes:    a 2D array containing the x,y and z positions of each node
%
% OUTPUT:
%   - spineNodes:       array containing the indexes of the spine nodes
%   - spineElements     array containing the indexes of the spine elements.
%                       In case there are multiple elements are sharing a 
%                       partof the spine, only the element with:
%                       - a third node having y=0 and z>0 AND
%                       - a fourth nodes having y>0
%                       is considered. 
%                       Such an element always exist and ist unique due to 
%                       the mesh partioning
%                       This allows to avoid applying the force mulitple
%                       times
%   - spineElementWeights:  array of length nElements, with 1 if element is
%                           part of spine, 0 else (according to the list
%                           spineElements from above).
%   - nodeIdxPosInElements: matrix of size nElements x 2. For the row
%                           corresponding to spine elements, it gives the
%                           column index (1 to 4) in the element of the 
%                           node close to the tail (1st column) and the 
%                           node close to the head (2nd column).
%
% Additional notes: -
%
% Last modified: 14/10/2023, Mathieu Dubied, ETH Zurich
function [spineNodes, spineElements, spineElementWeights, nodeIdxPosInElements] = find_spine_TET4(elements, nodes)
    nNodes = length(nodes(:,1));
    nElements = length(elements(:,1));

    % find spine nodes
    spineNodes = [];
    for n = 1:nNodes
        if nodes(n,2)==0 && nodes(n,3)==0
            spineNodes = [spineNodes; n];
        end
    end

    % find spine nodes in elements (logical matrix with 1 and 0)
    snInElements = ismember(elements,spineNodes);

    % select elements (row index) that have 2 spine nodes 
    sElCandidates = find(sum(snInElements,2) > 1);
    
    % select elements that have 2 spine nodes a third node with y=0, z>0
    spineElements = [];
    for c = 1:length(sElCandidates)
        % spine element candidate (elC) we check
        elC = sElCandidates(c);

        % get index of third and fourth nodes in candidate element  
        restNodesInElement = elements(elC,find(~snInElements(elC,:)));
        
        % get y position of third node in element
        if (nodes(restNodesInElement(1),2)==0 && ...
                nodes(restNodesInElement(1),3)>0 && ...
                nodes(restNodesInElement(2),2)>0) ...
            || ...
                (nodes(restNodesInElement(2),2)==0 && ...
                nodes(restNodesInElement(2),3)>0 && ...
                nodes(restNodesInElement(1),2)>0)
            spineElements = [spineElements; elC];
        end
    
    end

    % create logical vector for the spine elements (weight vector)
    spineElementWeights = zeros(nElements,1);
    spineElementWeights(spineElements) = 1;

    % create a matrix that gives the spine node column position in
    % "elements" for the elements that are spine elements
    nodeIdxPosInElements = zeros(nElements,2);
    for s = 1:length(spineElements)
        sEl = spineElements(s);                     % element number
        cIdx = find(snInElements(sEl,:));           % column index of spine nodes in this element
        sN = elements(sEl,snInElements(sEl,:));     % spine nodes numbers
        
        % check which nodes is closer to tail (x smaller)
        inv = 0;
        if nodes(sN(1),1) > nodes(sN(2),1)
            inv = 1;
        end

        % output column indexes [tail spine node, head spine node]
        if inv == 0
            snTail = cIdx(1);
            snHead = cIdx(2);
        else
            snTail = cIdx(2);
            snHead = cIdx(1);
        end

        nodeIdxPosInElements(sEl,:) = [snTail,snHead];        
    end
    

end