% gradient_cost_function_w_constraints_multi_objective
%
% Synthax:
% nablaLr = gradient_cost_function_w_constraints_multi_objective(xi,eta,s,V, ...
%                               AConstraint,bConstraint,barrierParam, wSize)
% Description:  gradient of the (reduced) cost function Lr for
%               multi-objective optimization
%
% Note:         Tailored to the case where the last parameter in the
%               parameter vector represent the actuation signal amplitude
%
% INPUTS: 
% (1) p:                current value for parameter vector p
% (2) eta:              solution for the reduced state variables
% (3) s:                solution for the sensitivity
% (4) V:                ROB matrix
% (5)-(6) A, b:         constraints on xi of the form Axi<b  
% (7) barrierParam:     parameter to scale (1/barrierParam) the barrier functions   
% (8) wSize:            windows size, cost is computed from 1-wSiue to N (included) 
% (9) w1:               weight for objective 1
% (10) w2:              weight for objective 2
% (11) alphaActu:       initial weight for objective 2
%                       
% OUTPUTS:
% (1) nablaLr:          gradient of the reduced cost function
%     
%
% Last modified: 07/02/2026, Mathieu Dubied, ETH Zürich

function nablaLr = gradient_cost_function_w_constraints_multi_objective(p,eta,s,V, ...
                               AConstraint,bConstraint,barrierParam, wSize,...
                               w1, w2, alphaActu)
    N = size(eta,2);
    nablaLr = zeros(size(p,1),1);
    nConstraints = size(bConstraint);
    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;
    
    for t=N-wSize:N
        % get gradient dfdxi_i (dfdp_i) a time step t        
        if size(p,1)>1
            % Objective 1
            s = double(s);
            dLdxi_i = -xDir.'*V*s(:,:,t);
            % Objective 2
            dLdpActu_i = alphaActu*[zeros(size(p(1:end-1)));1];
            
        else
            % Objective 1
            dLdxi_i = -xDir.'*V*s(:,t);
            
        end
        nablaLr = nablaLr + w1*dLdxi_i';
    end  
    
    % Objective 2
    dLdpActu = alphaActu*[zeros(size(p(1:end-1)));1];
    nablaLr = nablaLr + w2*dLdpActu;
    
    % part stemming from log barrier functions
    logBarrierD = zeros(size(p,1),1);
    for i = 1:nConstraints 
        logBarrierD = logBarrierD - 1/barrierParam(i)*1/(AConstraint(i,:)*p-bConstraint(i))*AConstraint(i,:).';
    end

    % final gradient
    nablaLr = nablaLr + logBarrierD;
  
end