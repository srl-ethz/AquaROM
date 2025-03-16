% reduced_cost_function_w_constraints_TET4
%
% Synthax:
% Lr = reduced_cost_function_w_constraints_TET4(N,tailProperties,spineProperties,eta0,eta,etad,etadd,xiRebuild,xi,dr,AConstraint,bConstraint,barrierParam)
%
% Description:  Computes the cost function value in the ROB, considering
%               constraints on the shape variation parameters xi. Upper and
%               lower bounds inequality are considered. These inequality
%               constraints are included as (log) barrier functions. 
%
% INPUTS: 
% (1) N:                number of time steps             
% (2) tailProperties:   properties of the tail pressure force
%                       (matrices, tail elements etc.)
% (3) spineProperties:  properties of the spine change in momentum
%                       (tensor, spine elements etc.)
% (4) x0:               initial node position in FOM
% (5) eta:              solution for the reduced state variables
% (6) etad:             solution for the reduced velocities
% (7) etadd:            solution for the reduced accelerations
% (8) xiRebuild:        current value for xi, after the last PROM rebuild
% (9) xi:               current value for xi, after first PROM build
% (10) dr:              reduced forward swimming direction vector
% (11)-(12) A, b:       constraints on xi of the form Axi<b  
% (12) barrierParam:    parameter to scale (1/barrierParam) the barrier functions      
%                   
%
% OUTPUTS:
% (1) L:   reduced cost function value
% (2) LwoB: cost function without barrier function
%     
%
% Last modified: 22/01/2023, Mathieu Dubied, ETH Zurich

function [L,LwoB] = reduced_cost_function_w_constraints_TET4(N,eta,xi,AConstraint,bConstraint,barrierParam,V)
    L = 0;
    LwoB = 0;
    nConstraints = size(bConstraint);
    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;
    
    for t=N-5:N
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