% reduced_cost_function_w_constraints_TET4
%
% Synthax:
% [L,LwoB] = reduced_cost_function_w_constraints_TET4(eta,xi,V,AConstraint,bConstraint,barrierParam,wSize)
 %   
% Description:  Computes the cost function value in the ROB, considering
%               constraints on the shape variation parameters xi. Upper and
%               lower bounds inequality are considered. These inequality
%               constraints are included as (log) barrier functions. 
%
% INPUTS: 
% (1) xi:               current value for xi, after first PROM build
% (2) eta:              solution for the reduced state variables
% (3) V:                ROB matrix
% (4)-(5) A, b:         constraints on xi of the form Axi<b  
% (6) barrierParam:     parameter to scale (1/barrierParam) the barrier functions   
% (7) wSize:            windows size, cost is computed from 1-wSiue to N (included)
%                   
% OUTPUTS:
% (1) L:   reduced cost function value
% (2) LwoB: cost function without barrier function
%     
% Last modified: 17/05/2025, Mathieu Dubied, ETH Zurich

function [L,LwoB] = reduced_cost_function_w_constraints_TET4(xi,eta,V,AConstraint,bConstraint,barrierParam,wSize)
    L = 0;
    LwoB = 0;
    N = size(eta,2);
    nConstraints = size(bConstraint);
    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;
    
    for t=N-wSize:N
        eta_i = eta(:,t);
        % add cost function at time step t to overall cost L       
        L = L - xDir.'*V*eta_i ;
    end
    
    % constraints (log barriers) to be included in the cost function
    logBarrier = 0;
    if nConstraints ~= 0
        for i = 1:nConstraints
            logBarrier = logBarrier - 1/barrierParam(i)*log(-AConstraint(i,:)*xi+bConstraint(i));
        end
    end
    
    % final cost
    LwoB = L;
    L = L + logBarrier;

    % print cost function without part stemming from barrier functions
    fprintf('Partial cost (w/o barrier): %.4f\n',LwoB)
    fprintf('Full cost (with barrier): %.4f\n',L)
end