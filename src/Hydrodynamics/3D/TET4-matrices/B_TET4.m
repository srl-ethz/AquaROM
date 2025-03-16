% return correct B as a function of the index of the tail node and the
% "head" node in the tail element. The "head" node is the node within the
% tail element that is closest to the head. The indexes are 1,2,3, or 4.
function B = B_TET4(tailNodeIdx,headNodeIdx)
    tailxIdx = tailNodeIdx*3-2;
    tailyIdx = tailNodeIdx*3-1;
    headxIdx = headNodeIdx*3-2;
    headyIdx = headNodeIdx*3-1;

    B = zeros(12,12); 
    B(tailxIdx,tailxIdx) = -1;
    B(tailyIdx,tailyIdx) = -1;
    B(tailxIdx,headxIdx) = 1;
    B(tailyIdx,headyIdx) = 1;
end