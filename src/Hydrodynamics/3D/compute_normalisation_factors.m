% compute_normalisation_factors
%
% Synthax:
% normalisationFactors = compute_normalisation_factors(nodes, elements, spineElements, nodeIdxPosInElements)
%
% Description: 
% returns normalisation factors, i.e., the inverse of the distance between
% each pair of spine nodes.
%
% INPUTS
%   - nodes:                a 2D array containing the x and y positions of 
%                           each node
%   - elements:             a 2D array in which each row contains the nodes' 
%                           indexes of an element
%   - spineElements:        array containing the indexes of the spine elements.
%                           In case there are two elements are sharing a face 
%                           on the spine, only the upper element is considered.
%                           This allows to avoid applying the force twice.
%   - nodeIdxPosInElements: matrix of size nElements x 2. For the row
%                           corresponding to spine elements, it give the
%                           column index (1 to 3) in the element of the 
%                           node close to the tail (1st column) and the 
%                           node close to the head (2nd column).
%
% OUTPUTS
%   - normalisationFactors: array of size nElements with zeros everywhere,
%                           except at the indexes of the spine elements.
%                           For these indexes, the inverse of the distance
%                           between the two spine nodes is returned.
%
% Additional notes: -
%
% Last modified: 06/10/2023, Mathieu Dubied, ETH Zurich

function normalisationFactors = compute_normalisation_factors(nodes, elements, spineElements, nodeIdxPosInElements)
    
    normalisationFactors = zeros(length(elements(:,1)),1);
    % loop over spine elements
    for s = 1:length(spineElements)
        sEl = spineElements(s);

        % get spine nodes of the element, one in tail direction, one in
        % head direction
        tailIdx = nodeIdxPosInElements(sEl,1);
        headIdx = nodeIdxPosInElements(sEl,2);
        tailPos = nodes(elements(sEl,tailIdx),:);
        headPos = nodes(elements(sEl,headIdx),:);

        % compute normalisation factor
        normalisationFactors(sEl) = 1/norm(headPos-tailPos);
    end
end
