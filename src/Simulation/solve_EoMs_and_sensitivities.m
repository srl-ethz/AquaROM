% solve_EoMs_and_sensitivities
%
% Synthax:
% TI_NL_PROM = solve_EoMs_and_sensitivities(V,ROM_Assembly,tailProperties,spineProperties,actuTop,actuBottom,h,tmax)
%
% Description: Computes the solutions for eta,dot{eta}, ddot{eta} for 
% [0,tmax] for a given ROM assembly and time step as well as the
% corresponding sensitivities (in a combined fashion)
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
% (8) kActu:                multiplicative factor for the actuation forces
% (9) h:                    time step for time integration
% (10) tmax:                simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_NL_PROM:           struct containing the solutions and related
%                           information
%     
% Last modified: 02/09/2024, Mathieu Dubied, ETH Zurich

function TI_NL_PROM = solve_EoMs_and_sensitivities(V,PROM_Assembly,fIntTensors,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax)

    fishDim = size(PROM_Assembly.Mesh.nodes,2);

    % SIMULATION PARAMETERS AND ICs _______________________________________
    eta0 = zeros(size(V,2),1);
    etad0 = zeros(size(V,2),1);
    etadd0 = zeros(size(V,2),1);

    U = PROM_Assembly.U;
    s0 = zeros(size(V,2),size(U,2));
    sd0 = zeros(size(V,2),size(U,2));
    sdd0 = zeros(size(V,2),size(U,2));

    % FORCES for EoMS: ACTUATION, REACTIVE FORCE, DRAG FORCE ______________
    % actuation force
    B1T = actuTop.B1;
    B1B = actuBottom.B1;
    B2T = actuTop.B2;
    B2B = actuBottom.B2;
    
    actuSignalT = @(t) kActu/2*(-0.2*sin(t*2*pi));    % to change below as well if needed
    actuSignalB = @(t) kActu/2*(0.2*sin(t*2*pi));
    
    fActu = @(t,q)  kActu/2*(-0.2*sin(t*2*pi))*(B1T+B2T*q) + ...
                    kActu/2*(0.2*sin(t*2*pi))*(B1B+B2B*q);

    % tail pressure force properties
    A = tailProperties.A;
    B = tailProperties.B;
    R = tailProperties.R;
    wTail = tailProperties.w;
    VTail = tailProperties.V;
    UTail = tailProperties.U;
    nodes = PROM_Assembly.Mesh.nodes;
    iDOFs = tailProperties.iDOFs;

    if fishDim == 3
        tailProperties.mTilde = 0.25*pi*1000*(tailProperties.z*2)^2;
    end

    x0 = reshape(PROM_Assembly.Mesh.nodes.',[],1);  
                       
    fTail = @(q,qd)  0.5*tailProperties.mTilde*wTail^3*VTail.'*(dot(A*VTail*qd,R*B*(x0(tailProperties.iDOFs)+VTail*q))).^2* ...
                        B*(x0(tailProperties.iDOFs)+VTail*q);
    
    % spine change in momentum (solve for xi=0, so that we do not need f1
    % and f2 terms, nor the terms in U)
    if fishDim == 3
        Txx2 = spineProperties.tensors.Txx.f0;
        TxV3 = spineProperties.tensors.TxV.f0;
        TVx3 = spineProperties.tensors.TVx.f0;
        TVV4 = spineProperties.tensors.TVV.f0;
        
        fSpine = @(q,qd,qdd) double(Txx2)*qdd ...
            + double(ttv(ttv(TxV3,q,3),qdd,2) ...
            + ttv(ttv(TVx3,q,3),qdd,2) ...
            + ttv(ttv(ttv(TVV4,q,4),q,3),qdd,2) ...
            + ttv(ttv(TVx3,qd,3),qd,2) ...
            + ttv(ttv(ttv(TVV4,q,4),qd,3),qd,2) ...
            + ttv(ttv(TxV3,qd,3),qd,2) ...
            + ttv(ttv(ttv(TVV4,qd,4),q,3),qd,2));
    else
        Txx = spineProperties.tensors.Txx;
        TxV = spineProperties.tensors.TxV;
        TVx = spineProperties.tensors.TVx;
        TVV = spineProperties.tensors.TVV;
        
        fSpine = @(q,qd,qdd) double(Txx)*qdd ...
            + double(ttv(ttv(TxV,q,3),qdd,2) ...
            + ttv(ttv(TVx,q,3),qdd,2) ...
            + ttv(ttv(ttv(TVV,q,4),q,3),qdd,2) ...
            + ttv(ttv(TVx,qd,3),qd,2) ...
            + ttv(ttv(ttv(TVV,q,4),qd,3),qd,2) ...
            + ttv(ttv(TxV,qd,3),qd,2) ...
            + ttv(ttv(ttv(TVV,qd,4),q,3),qd,2));
    end

    % drag force
    T3 = dragProperties.tensors.Tr3;
    fDrag = @(qd) double(ttv(ttv(T3,qd,3),qd,2));

    % FORCES DERIVATIVES FOR SENSITIVITY ANALYSIS _________________________
    % internal forces
    pd_fint = @(eta)DpROM_derivatives_lightweight(eta,fIntTensors);

    % actuation forces
    pd_actuTop = @(a)PROM_actu_derivatives(actuTop,a);
    pd_actuBottom = @(a)PROM_actu_derivatives(actuBottom,a);

    % tail pressure force
    xi=zeros(size(PROM_Assembly.U,2),1);

    if fishDim == 3
        z0 = tailProperties.z;
        Uz = tailProperties.Uz;
        pd_tail = @(eta,etad) PROM_tail_pressure_derivatives_TET4_lightweight(eta,etad,A,B,R,wTail,x0(iDOFs),VTail,UTail,z0,Uz);
    else 
        mTilde = tailProperties.mTilde;
        pd_tail = @(eta,etad) PROM_tail_pressure_derivatives(eta,etad,A,B,R,mTilde,wTail,x0(iDOFs),xi,VTail,UTail);                                                           
    end

    % spine change in momentum
    spineTensors = spineProperties.tensors;
    if fishDim == 3
        pd_spine = @(eta,etad,etadd)PROM_spine_momentum_derivatives_TET4_lightweight(eta,etad,etadd,spineTensors);
    else
        pd_spine = @(eta,etad,etadd)PROM_spine_momentum_derivatives(eta,etad,etadd,xi,spineTensors);
    end

    % drag force
    dragTensors = dragProperties.tensors;
    pd_drag = @(etad) PROM_drag_derivatives_lightweight(etad,dragTensors);


    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object for nonlinear time integration
    TI_NL_PROM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-6,'combinedSensitivity',true);
    
    % modal nonlinear Residual evaluation function handle
    Residual_NL_red = @(eta,etad,etadd,t)residual_reduced_nonlinear_actu_hydro_PROM(eta,etad,etadd, ...
        t,PROM_Assembly,fIntTensors,fActu,fTail,fSpine,fDrag,actuTop,actuBottom,actuSignalT,actuSignalB,tailProperties,spineProperties,dragProperties,R,x0);

    % residual function handle
    Residual_sens = @(s,sd,sdd,etaSol,etadSol,etaddSol,drdqdd, drdqd, drdq,t)residual_linear_sens_combined(s,sd,sdd,t, ...
        etaSol,etadSol,etaddSol, drdqdd, drdqd, drdq,...
        pd_fint,pd_tail,pd_spine,pd_drag,pd_actuTop,pd_actuBottom, ...
        actuSignalT,actuSignalB);

    % time integration 
    TI_NL_PROM.Integrate(eta0,etad0,etadd0,tmax,Residual_NL_red, ...
        'ResidualSens',Residual_sens,'s0',s0,'sd0',sd0,'sdd0',sdd0);
    TI_NL_PROM.Solution.u = V * TI_NL_PROM.Solution.q; % get full order solution

end 
