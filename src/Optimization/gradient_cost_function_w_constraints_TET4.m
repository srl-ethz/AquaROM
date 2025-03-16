% gradient_cost_function_w_constraints_TET4
%
% Synthax:
% nabla_Lr = gradient_cost_function_w_constraints_TET4(dr,xiRebuild,xi,eta0,eta,etad,etadd,s,sd,sdd,tailProperties,spineProperties,AConstraint,bConstraint,barrierParam)
%
% Description:  gradient of the (reduced) cost function Lr. The gradient is
%               analytical and based on the hydrodynamic tensors
%
% INPUTS: 
% (1) dr:               reduced forward swimming direction vector
% (2) xiRebuild:        current value for xi, after the last PROM rebuild
% (3) xi:               current value for xi, after first PROM build
% (4) x0:               initial node position in FOM
% (5) eta:              solution for the reduced state variables
% (6) etad:             solution for the reduced velocities
% (7) etadd:            solution for the reduced accelerations
% (8) s:                solution for the sensitivity
% (9) sd:               solution for the sensitivity derivative
% (10) sdd:             solution for the sensitivity 2nd derivative
% (11)tailProperties:   properties of the tail pressure force
%                       (matrices, tail elements etc.)
% (12) spineProperties: properties of the spine change in momentum
%                       (tensor, spine elements etc.)
% (13)-(14) A,b:        constraints on xi of the form Axi<b 
% (14) barrierParam:    parameter to scale (1/barrierParam) the barrier functions 
%                       
% OUTPUTS:
% (1) nablaLr:          gradient of the reduced cost function
%     
%
% Last modified: 22/01/2023, Mathieu Dubied, ETH Zürich

function nablaLr = gradient_cost_function_w_constraints_TET4(xi,eta,s,AConstraint,bConstraint,barrierParam,V)
    N = size(eta,2);
    nablaLr = zeros(size(xi,1),1);
    nConstraints = size(bConstraint);
    xDir = zeros(size(V,1),1);
    xDir(1:3:end) = 1;
    
    for t=N-5:N
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