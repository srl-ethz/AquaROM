% compute_AB_tail_pressure
%
% Synthax:
% force = compute_AB_tail_pressure(nodeIdxPosInTailElement)
%
% Description: 
% Computes the A and B matrices needed to express the tail pressure force.
%
% INPUTS
%   - nodeIdxPosInTailElement: array of size 1 x 2. It gives the
%                               column index (1 to 3) in the element of the 
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
% Last modified: 10/10/2023, Mathieu Dubied, ETH Zurich
function [A,B] = compute_AB_tail_pressure(nodeIdxPosInTailElement)
    % get matrix corresponding to the configuration
    if nodeIdxPosInTailElement(1)==1         % conf 1-2 and 1-3
        A = A_conf1;
        if nodeIdxPosInTailElement(2)==2
            B = B_conf1_2;
        else
            B = B_conf1_3;
        end
    elseif nodeIdxPosInTailElement(1)==2     % conf 2-1 and 2-3
        A = A_conf2;
        if nodeIdxPosInTailElement(2)==1
            B = B_conf2_1;
        else
            B = B_conf2_3;
        end
    else                                    % conf 3-1 and 3-2
        A = A_conf3;
        if nodeIdxPosInTailElement(2)==1
            B = B_conf3_1;
        else
            B = B_conf3_2;
        end
    end

end