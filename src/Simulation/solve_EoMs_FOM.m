% solve_EoMs_FOM
%
% Synthax:
% TI_NL_FOM = solve_EoMs_FOM(Assembly, elements, tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax)
%
% Description: Computes the solutions for u, \dot{u}, \ddot{u} for 
% [0,tmax] for a given FOM assembly and time step
%
% INPUTS:  
% (1) Assembly:             FOM assembly    
% (2) elements              elements of the assembly
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
% (1) TI_NL_FOM:            struct containing the solutions and related
%                           information
%     
% Last modified: 17/10/2024, Mathieu Dubied, ETH Zurich

function TI_NL_FOM = solve_EoMs_FOM(Assembly, elements, tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax)
    
    % SIMULATION PARAMETERS AND ICs _______________________________________
    nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);
    nDOFs = Assembly.Mesh.nDOFs;
    q0 = zeros(nUncDOFs,1);
    qd0 = zeros(nUncDOFs,1);
    qdd0 = zeros(nUncDOFs,1);

    % FORCES: ACTUATION, REACTIVE FORCE, DRAG FORCE _______________________
    % actuation force
    B1T = actuTop.B1;
    B1B = actuBottom.B1;
    B2T = actuTop.B2;
    B2B = actuBottom.B2; 
    
    actuSignalT = @(t) kActu/2*(-0.2*sin(t*2*pi));    % to change below as well if needed
    actuSignalB = @(t) kActu/2*(0.2*sin(t*2*pi));
    
    fActu = @(t,q)  kActu/2*(-0.2*sin(t*2*pi))*(B1T+B2T*q) + ...
                    kActu/2*(0.2*sin(t*2*pi))*(B1B+B2B*q);
%                 
%     fActu = @(t,q) kActu*tip_actuation_force_FOM(t,q,Assembly,elements,tailProperties);
%     actuSignalT = kActu;
%     fActu = @(t,q) kActu*tip_actuation_force_FOM(t,q,Assembly,elements,tailProperties);
    
%     fActu = @(t) tip_actuation_force_FOM_1(t,kActu,Assembly,tailProperties);
    
    % tail pressure force
    fTail = @(q,qd) 1*tail_force_FOM(q,qd,Assembly,elements,tailProperties);

    % spine force
    fSpine = @(q,qd,qdd) 1*spine_force_FOM(q,qd,qdd,Assembly,elements,spineProperties);

    % drag force
    fDrag = @(qd) 1*drag_force_FOM(qd,dragProperties);

    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate object for nonlinear time integration
    TI_NL_FOM = ImplicitNewmark('timestep',h,'alpha',0.005,'MaxNRit',60,'RelTol',1e-5);

    % modal nonlinear Residual evaluation function handle
    Residual_NL = @(q,qd,qdd,t)residual_nonlinear_actu_hydro(q,qd,qdd, ...
        t,Assembly,fActu,fTail,fSpine,fDrag,actuTop,actuBottom,actuSignalT,actuSignalB,tailProperties,spineProperties,elements);

    % time integration 
    TI_NL_FOM.Integrate(q0,qd0,qdd0,tmax,Residual_NL);
    TI_NL_FOM.Solution.u = Assembly.unconstrain_vector(TI_NL_FOM.Solution.q);

end 
