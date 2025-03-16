% reduced_constant_vector
%
% Synthax:
% dr = reduced_constant_vector(d,V)
%
% Description: Convert a constant direction vector to its reduced version
%
% INPUTS: 
% (1) d:    2D forward swimming direction vector             
% (2) V:    ROB
%
%
% OUTPUTS:
% (1) dr:   reduced order version of d    
%     
%
% Additional notes:
%   - specific implementation for TRI elements
%
% Last modified: 24/10s/2023, Mathieu Dubied, ETH Zürich

function dr = reduced_constant_vector(d,V,dim)
    n = size(V,1);
    m = size(V,2);
    dr = zeros(m,1);
%     for i=1:6:n
%         Ve = V(i:i+5,:);
%         de = [d;d;d];
%         dr = dr + Ve.'*de;  
%     end
    if dim == 3
        dFull = repmat(d',1,n/3)';
    else
        dFull = repmat(d',1,n/2);
    end
    dr = V'*dFull;
end
