% reduced_cost_function_w_constraints_multi_objectives
%
% Synthax:
% [L,LwoB,LObj1,LObj2] = reduced_cost_function_w_constraints_TET4(eta,xi,V,AConstraint,bConstraint,barrierParam,wSize)
 %   
% Description:  Computes the cost function value in the ROB, considering
%               constraints on the shape variation parameters xi. Upper and
%               lower bounds inequality are considered. These inequality
%               constraints are included as (log) barrier functions.
%
% Note:         Tailored to the case where the last parameter in the
%               parameter vector represent the actuation signal amplitude
%
% INPUTS: 
% (1) p:                current value for parameter vector p
% (2) eta:              solution for the reduced state variables
% (3) V:                ROB matrix
% (4)-(5) A, b:         constraints on xi of the form Axi<b  
% (6) barrierParam:     parameter to scale (1/barrierParam) the barrier functions   
% (7) wSize:            windows size, cost is computed from 1-wSiue to N (included)
% (8) w1:               weight for objective 1
% (9) w2:               weight for objective 2
% (10) alphaActu:       initial weight for objective 2
%                   
% OUTPUTS:
% (1) L:        reduced cost function value
% (2) LwoB:     cost function without barrier function
% (3) LObj1:    cost with respect to objective 1
% (4) LObj2:    cost with respect to objective 2
%     
% Last modified: 07/02/2026, Mathieu Dubied, ETH Zurich

function [L,LwoB,LObj1,LObj2] = ...
    reduced_cost_function_w_constraints_multi_objective(p,eta,V,...
                        AConstraint,bConstraint,barrierParam,wSize, ...
                        w1, w2, alphaActu)
    L = 0;
    LObj1 = 0;
    LObj2 = 0;
    LwoB = 0;
    N = size(eta,2);
    nConstraints = size(bConstraint);
    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;
    
    for t=N-wSize:N
        eta_i = eta(:,t);
        % Objective 1: maximize swimming distance     
        LObj1 = LObj1 - xDir.'*V*eta_i ;
    end
    
    % Objective 2: minimize energy use, i.e., actuation amplitude
    LObj2 = LObj2 + alphaActu*p(end);
    
    L = w1*LObj1 + w2*LObj2;
    
    % constraints (log barriers) to be included in the cost function
    logBarrier = 0;
    if nConstraints ~= 0
        for i = 1:nConstraints
            logBarrier = logBarrier - 1/barrierParam(i)*log(-AConstraint(i,:)*p+bConstraint(i));
        end
    end
    
    % final cost
    LwoB = L;
    L = L + logBarrier;

    % print cost function without part stemming from barrier functions
    fprintf('Objective 1 cost: %.4f\n',LObj1)
    fprintf('Objective 2 cost: %.4f\n',LObj2)
    fprintf('Objective 1 weighted cost: %.4f\n',w1*LObj1)
    fprintf('Objective 2 weighted cost: %.4f\n',w2*LObj2)
    fprintf('Partial cost (w/o barrier): %.4f\n',LwoB)
    fprintf('Full cost (with barrier): %.4f\n',L)
end