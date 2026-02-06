% gradient_cost_function_wo_constraints_TET4
%
% Synthax:
% nablaLr = gradient_cost_function_wo_constraints_TET4(xi,eta,s,V, ...
%                               wSize)
% Description:  gradient of the (reduced) cost function Lr. The gradient is
%               analytical and based on the hydrodynamic tensors
%
% INPUTS: 
% (1) xi:               current value for xi, after first PROM build
% (2) eta:              solution for the reduced state variables
% (3) s:                solution for the sensitivity
% (4) V:                ROB matrix   
% (8) wSize:            windows size, cost is computed from 1-wSiue to N (included) 
%                       
% OUTPUTS:
% (1) nablaLr:          gradient of the reduced cost function
%     
%
% Last modified: 05/02/2025, Mathieu Dubied, ETH Zürich

function nablaLr = gradient_cost_function_wo_constraints_TET4(xi,eta,s,V,wSize)
    N = size(eta,2);
    nablaLr = zeros(size(xi,1),1);
    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;
    
    for t=N-wSize:N
        % get gradient dfdxi_i (dfdp_i) a time step t        
        if size(xi,1)>1
            s = double(s);
            dLdxi_i = -xDir.'*V*s(:,:,t);
        else
            dLdxi_i = -xDir.'*V*s(:,t);
        end
        
        % add to overall gradient
        nablaLr = nablaLr + dLdxi_i';
    end  
  
end