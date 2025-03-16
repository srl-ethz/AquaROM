% solve_EoMs_and_sensitivities_actu
%
% Synthax:
% TI_NL_PROM = solve_EoMs_and_sensitivities_actu(V,ROM_Assembly,tailProperties,spineProperties,actuTop,actuBottom,h,tmax)
%
% Description: Computes the solutions for eta,dot{eta}, ddot{eta} for 
% [0,tmax] for a given ROM assembly and time step as well as the
% corresponding sensitivities of the actuation (in a combined fashion)
%
% INPUTS: 
% (1) V:                    ROB   
% (2) ROM_Assembly:         ROM assembly           
% (3) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (4) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (5) dragProperties:       properties of the drag force
% (6) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (7) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
% (8) h:                    time step for time integration
% (9) tmax:                 simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_NL_PROM:           struct containing the solutions and related
%                           information
%     
% Last modified: 10/11/2023, Mathieu Dubied, ETH Zurich

function TI_NL_PROM = solve_EoMs_and_sensitivities_actu(V,PROM_Assembly,fIntTensors,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,p)

    fishDim = size(PROM_Assembly.Mesh.nodes,2);

    % SIMULATION PARAMETERS AND ICs _______________________________________
    eta0 = zeros(size(V,2),1);
    etad0 = zeros(size(V,2),1);
    etadd0 = zeros(size(V,2),1);

    s0 = zeros(size(V,2),length(p));
    sd0 = zeros(size(V,2),length(p));
    sdd0 = zeros(size(V,2),length(p));

    % FORCES for EoMS: ACTUATION, REACTIVE FORCE, DRAG FORCE ______________
    % actuation force
    B1T = actuTop.B1;
    B1B = actuBottom.B1;
    B2T = actuTop.B2;
    B2B = actuBottom.B2;
    k = 3.0; % k=400 for 3D
    
    actuSignalT = @(t) actuation_signal_1(k,t,p);    
    actuSignalB = @(t) -actuation_signal_1(k,t,p);
    
    fActu = @(t,q) actuation_force_1(k,t,q,B1T,B2T, B1B,B2B,p);

    % tail pressure force properties
    A = tailProperties.A;
    B = tailProperties.B;
    R = tailProperties.R;
    wTail = tailProperties.w;
    VTail = tailProperties.V;

    if fishDim == 3
        tailProperties.mTilde = 0.25*pi*1000*(tailProperties.z*2)^2;
    end

    x0 = reshape(PROM_Assembly.Mesh.nodes.',[],1);  
                       
    fTail = @(q,qd)  0.5*tailProperties.mTilde*wTail^3*VTail.'*(dot(A*VTail*qd,R*B*(x0(tailProperties.iDOFs)+VTail*q))).^2* ...
                        B*(x0(tailProperties.iDOFs)+VTail*q);
    
    % spine change  momentum
%     Txx = spineProperties.tensors.Txx;
%     TxV = spineProperties.tensors.TxV;
%     TVx = spineProperties.tensors.TVx;
%     TVV = spineProperties.tensors.TVV;
%     
%     fSpine = @(q,qd,qdd) double(Txx)*qdd ...
%         + double(ttv(ttv(TxV,q,3),qdd,2) ...
%         + ttv(ttv(TVx,q,3),qdd,2) ...
%         + ttv(ttv(ttv(TVV,q,4),q,3),qdd,2) ...
%         + ttv(ttv(TVx,qd,3),qd,2) ...
%         + ttv(ttv(ttv(TVV,q,4),qd,3),qd,2) ...
%         + ttv(ttv(TxV,qd,3),qd,2) ...
%         + ttv(ttv(ttv(TVV,qd,4),q,3),qd,2));
    if isfield(spineProperties.tensors,'f0') == 1   % PROM was constructed
        Txx2 = spineProperties.tensors.Txx.f0;
        TxV3 = spineProperties.tensors.TxV.f0;
        TVx3 = spineProperties.tensors.TVx.f0;
        TVV4 = spineProperties.tensors.TVV.f0;
    else                                            % ROM was constructed
        Txx2 = spineProperties.tensors.Txx; 
        TxV3 = spineProperties.tensors.TxV;
        TVx3 = spineProperties.tensors.TVx;
        TVV4 = spineProperties.tensors.TVV;
    end

    fSpine = @(q,qd,qdd) double(Txx2)*qdd ...
        + double(ttv(ttv(TxV3,q,3),qdd,2) ...
        + ttv(ttv(TVx3,q,3),qdd,2) ...
        + ttv(ttv(ttv(TVV4,q,4),q,3),qdd,2) ...
        + ttv(ttv(TVx3,qd,3),qd,2) ...
        + ttv(ttv(ttv(TVV4,q,4),qd,3),qd,2) ...
        + ttv(ttv(TxV3,qd,3),qd,2) ...
        + ttv(ttv(ttv(TVV4,qd,4),q,3),qd,2));


    % drag force
    T3 = dragProperties.tensors.Tr3;
    fDrag = @(qd) double(ttv(ttv(T3,qd,3),qd,2));

    % ACTUATION FORCES DERIVATIVES FOR SENSITIVITY ANALYSIS _______________
    % actuation forces
    pd_actu = @(t,q,p)derivatives_actuation_force_1(k,t,q,B1T,B2T,B1B,B2B,p);
        
    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object for nonlinear time integration
    TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6,'combinedSensitivity',true);
    
    % modal nonlinear Residual evaluation function handle
    Residual_NL_red = @(eta,etad,etadd,t)residual_reduced_nonlinear_actu_hydro(eta,etad,etadd, ...
        t,PROM_Assembly,fIntTensors,fActu,fTail,fSpine,fDrag,actuTop,actuBottom,actuSignalT,actuSignalB,tailProperties,spineProperties,dragProperties,R,x0);

    % residual function handle
    Residual_sens = @(s,sd,sdd,t,etaSol,drdqdd, drdqd, drdq)residual_linear_sens_combined_actu(s,sd,sdd,t,etaSol,drdqdd,drdqd,drdq,pd_actu,p);                                                                                     
    
    % time integration 
    TI_NL_PROM.Integrate(eta0,etad0,etadd0,tmax,Residual_NL_red, ...
        'ResidualSens',Residual_sens,'s0',s0,'sd0',sd0,'sdd0',sdd0,'actuOnly',true);
    TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

end 
