% gradient_cost_function_w_constraints
%
% Synthax:
% nabla_Lr = gradient_cost_function_w_constraints(dr,xiRebuild,xi,eta0,eta,etad,etadd,s,sd,sdd,tailProperties,spineProperties,AConstraint,bConstraint,barrierParam)
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
% Last modified: 15/10/2023, Mathieu Dubied, ETH Zürich

function nablaLr = gradient_cost_function_w_constraints(dr,xiRebuild,xi,x0,eta,etad,etadd,s,sd,sdd,tailProperties,spineProperties,AConstraint,bConstraint,barrierParam,V)
    N = size(eta,2);
    nablaLr = zeros(size(xi,1),1);
    nConstraints = size(bConstraint);
 
    %barrierParam = 14000;%400; for C:400

    % tail pressure force properties
    A = tailProperties.A;
    B = tailProperties.B;
    R = [0 -1 0 0 0 0;
         1 0 0 0 0 0;
         0 0 0 -1 0 0;
         0 0 1 0 0 0;
         0 0 0 0 0 1;
         0 0 0 0 -1 0];     % 90 degrees rotation counterclock-wise
    wTail = tailProperties.w;
    VTail = tailProperties.V;
    UTail = tailProperties.U;
    mTilde = tailProperties.mTilde;
    x0Tail = x0(tailProperties.iDOFs);
   
    % spine change in momentum
    tensorsSpine = spineProperties.tensors;


    xDir = zeros(size(V,1),1);
    xDir(1:2:end) = 1;
    
    for t=1:N-2 
        % tail pressure force
        % derTail = PROM_tail_pressure_derivatives(eta(:,t),etad(:,t),A,B,R,mTilde,wTail,x0Tail,xiRebuild,VTail,UTail); 
        % dfTaildq = derTail.dfdq;               
        % dfTaildqd = derTail.dfdqd;
        % dfTaildp = derTail.dfdp;
        % 
        % % spine change in momentum
        % derSpine = PROM_spine_momentum_derivatives(eta(:,t),etad(:,t),etadd(:,t),xiRebuild,tensorsSpine);        
        % dfSpinedq = derSpine.dfdq;  
        % dfSpinedqd = derSpine.dfdqd;
        % dfSpinedqdd = derSpine.dfdqdd;
        % dfSpinedp = derSpine.dfdp;

        % get gradient dfdxi_i (dfdp_i)         
        % if size(xi,1)>1
        %     s = double(s);
        %     sd = double(sd);
        %     sdd = double(sdd);
        %     dfdxi_i = dfTaildp + dfTaildq*s(:,:,t) + dfTaildqd*sd(:,:,t) ...
        %             + dfSpinedp + dfSpinedq*s(:,:,t) ...
        %             + dfSpinedqd*sd(:,:,t) + dfSpinedqdd*sdd(:,:,t); 
        %     % dLdxi_i = -[1;0;1;0;1;0].'*VTail*s(:,:,t);
        % 
        % else
        %     dfdxi_i = dfTaildp + dfTaildq*s(:,t) + dfTaildqd*sd(:,t) ...
        %             + dfSpinedp + dfSpinedq*s(:,t) ...
        %             + dfSpinedqd*sd(:,t) + dfSpinedqdd*sdd(:,t); 
        % 
        %     % dLdxi_i = -[1;0;1;0;1;0].'*VTail*s(:,t);
        % end

        if size(xi,1)>1
            s = double(s);
            % sd = double(sd);
            % sdd = double(sdd);
            % dfdxi_i = dfTaildp + dfTaildq*s(:,:,t) + dfTaildqd*sd(:,:,t) ...
            %         + dfSpinedp + dfSpinedq*s(:,:,t) ...
            %         + dfSpinedqd*sd(:,:,t) + dfSpinedqdd*sdd(:,:,t); 
            dLdxi_i = -xDir.'*V*s(:,:,t);

        else
            % dfdxi_i = dfTaildp + dfTaildq*s(:,t) + dfTaildqd*sd(:,t) ...
            %         + dfSpinedp + dfSpinedq*s(:,t) ...
            %         + dfSpinedqd*sd(:,t) + dfSpinedqdd*sdd(:,t); 

            dLdxi_i = -xDir.'*V*s(:,t);
        end
        
        
        % part stemming from log barrier functions
        logBarrierDInTimeStep = zeros(size(xi,1),1);
       
        for i = 1:nConstraints 
            logBarrierDInTimeStep = logBarrierDInTimeStep - 1/barrierParam*1/(AConstraint(i,:)*xi-bConstraint(i))*AConstraint(i,:).';
        end
        

        % final gradient
        % nablaLr = nablaLr - (dr.'*dfdxi_i).' + logBarrierDInTimeStep;
        % nablaLr = nablaLr - dLdxi_i + logBarrierDInTimeStep;
        nablaLr = nablaLr + dLdxi_i' + logBarrierDInTimeStep;
    end  
  
end