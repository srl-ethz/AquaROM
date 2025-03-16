% solve_EoMs_FOM_RK4
%
% Synthax:
% TI_NL_FOM = solve_EoMs_FOM_RK4(Assembly, elements, tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax)
%
% Description: Computes the solutions for u, \dot{u}, \ddot{u} for 
% [0,tmax] for a given FOM assembly and time step, using RK4
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
% Last modified: 11/11/2024, Mathieu Dubied, ETH Zurich

function TI_NL_FOM = solve_EoMs_FOM_RK4(Assembly, elements, tailProperties,spineProperties,dragProperties,actuTop,actuBottom,kActu,h,tmax)
    
    % SIMULATION PARAMETERS AND ICs _______________________________________
    nUncDOFs = size(Assembly.Mesh.EBC.unconstrainedDOFs,2);
    q0 = zeros(nUncDOFs,1);     % we solve for the displacement, not the position
    qd0 = zeros(nUncDOFs,1);
    y0 = [q0;qd0];

    % FORCES: ACTUATION, REACTIVE FORCE, DRAG FORCE _______________________
    % actuation force 
    B1T = actuTop.B1;
    B1B = actuBottom.B1;
    B2T = actuTop.B2;
    B2B = actuBottom.B2; 
        
    fActu = @(t,q)  kActu/2*(-0.2*sin(t*2*pi))*(B1T+B2T*q) + ...
                    kActu/2*(0.2*sin(t*2*pi))*(B1B+B2B*q);
    
    % tail pressure force
    fTail = @(q,qd) 0*tail_force_FOM(q,qd,Assembly,elements,tailProperties);

    % spine force
    fSpine = @(q,qd,qdd) 0*spine_force_FOM(q,qd,qdd,Assembly,elements,spineProperties);

    % drag force
    fDrag = @(qd) 0*drag_force_FOM(qd,dragProperties);

    % NONLINEAR TIME INTEGRATION __________________________________________
    
    % instantiate Runge-Kutta solver
    TI_NL_FOM = RungeKutta4('timestep', h);
    
    % define dynamics and convert it to a first order ODE
    Dynamics_RK4 = @(y,t) dynamics_for_RK4(y,t,Assembly,fActu,fTail,fSpine,fDrag);
    
    % time integration
    TI_NL_FOM.Integrate(y0, tmax, Dynamics_RK4);
    
    % store results/convert to nodes' position
    TI_NL_FOM.Solution.u = Assembly.unconstrain_vector(TI_NL_FOM.Solution.y(1:length(TI_NL_FOM.Solution.y)/2,:));


end 
