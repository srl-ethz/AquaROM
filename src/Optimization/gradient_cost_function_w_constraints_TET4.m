% gradient_cost_function_w_constraints_TET4
%
% Synthax:
% nablaLr = gradient_cost_function_w_constraints_TET4(xi,eta,s,V, ...
%                               AConstraint,bConstraint,barrierParam, wSize)
% Description:  gradient of the (reduced) cost function Lr. The gradient is
%               analytical and based on the hydrodynamic tensors
%
% INPUTS: 
% (1) xi:               current value for xi, after first PROM build
% (2) eta:              solution for the reduced state variables
% (3) s:                solution for the sensitivity
% (4) V:                ROB matrix
% (5)-(6) A, b:         constraints on xi of the form Axi<b  
% (7) barrierParam:     parameter to scale (1/barrierParam) the barrier functions   
% (8) wSize:            windows size, cost is computed from 1-wSiue to N (included) 
%                       
% OUTPUTS:
% (1) nablaLr:          gradient of the reduced cost function
%     
%
% Last modified: 17/05/2025, Mathieu Dubied, ETH Zürich

function nablaLr = gradient_cost_function_w_constraints_TET4(xi,eta,s,V, ...
                               AConstraint,bConstraint,barrierParam, wSize)
    N = size(eta,2);
    nablaLr = zeros(size(xi,1),1);
    nConstraints = size(bConstraint);
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
    
    % part stemming from log barrier functions
    logBarrierD = zeros(size(xi,1),1);
    for i = 1:nConstraints 
        logBarrierD = logBarrierD - 1/barrierParam(i)*1/(AConstraint(i,:)*xi-bConstraint(i))*AConstraint(i,:).';
    end

    % final gradient
    nablaLr = nablaLr + logBarrierD;
  
end