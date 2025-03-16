% return correct A as a function of the index of the tail node in the tail 
% element. The indexes are 1,2,3, or 4.
function A = A_TET4(tailNodeIdx)
    tailxIdx = tailNodeIdx*3-2;
    tailyIdx = tailNodeIdx*3-1;

    A = zeros(12,12); 
    A(tailxIdx,tailxIdx) = 1;
    A(tailyIdx,tailyIdx) = 1;
end