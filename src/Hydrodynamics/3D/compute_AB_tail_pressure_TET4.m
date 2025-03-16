% compute_AB_tail_pressure_3D
%
% Synthax:
% force = compute_AB_tail_pressure_TET4(nodeIdxPosInTailElement)
%
% Description: 
% Computes the A and B matrices needed to express the tail pressure force.
%
% INPUTS
%   - nodeIdxPosInTailElement:  array of size 1 x 2. It gives the
%                               column index (1 to 4) in the element of the 
%                               node close to the tail (1st column) and the 
%                               node close to the head (2nd column).

% OUTPUTS
%   - A:                        matrix needed to express the tail pressure
%                               force
%   - B:                        matrix needed to express the tail pressure
%                               force
%       
%
% Additional notes: -
%
% Last modified: 16/10/2023, Mathieu Dubied, ETH Zurich
function [A,B] = compute_AB_tail_pressure_TET4(nodeIdxPosInTailElement)

    % get matrix corresponding to the configuration
    A = A_TET4(nodeIdxPosInTailElement(1));
    B = B_TET4(nodeIdxPosInTailElement(1),nodeIdxPosInTailElement(2));

end







