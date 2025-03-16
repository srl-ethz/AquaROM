% find_fixed_point
%
% Synthax:
% fixedPointPos = find_fixed_point(nodes, desiredFixedPoint)
%
% Description: 
% Find the point of the spine that is the closest to the fixed point
%
% INPUTS
%   - 
%
% OUTPUT:
% Additional notes: -
%
% Last modified: 03/02/2024, Mathieu Dubied, ETH Zurich
function fixedPointPos = find_fixed_point(nodes, desiredFixedPoint)
    
    % select nodes on the spine
    nodeSelection = [];
    for n = 1:size(nodes,1)
        if nodes(n,2) == 0 && nodes(n,3) == 0
            nodeSelection = [nodeSelection,n];
        end
    end
   
    % find index in nodeSelection which correspond to the node closest to
    % desiredFixedPoint
    dist = abs(nodes(nodeSelection,1) - desiredFixedPoint);
    [~,ind]=min(dist);
    
    % find corresponding node x-position
    fixedPointPos = nodes(nodeSelection(ind),1);
end