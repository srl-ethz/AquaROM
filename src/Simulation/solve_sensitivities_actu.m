% solve_sensitivities_actu
%
% Synthax:
% TI_sens = solve_sensitivities_actu(V,xi_k,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,etaSol,etadSol,etaddSol,h,tmax)
%
% Description: Computes the solutions for the sensitivity S for [0,tmax] 
% for a given ROM assembly, actuation parameter p_k, and time step h
%
% INPUTS: 
% (1) V:                    ROB 
% (2) xi_k:                 current shape variation xi to consider
% (3) PROM_Assembly:        PROM assembly  
% (4) tensors_PROM          internal forces tensors
% (5) tailProperties:       properties of the tail pressure force
%                           (matrices, tail elements etc.)
% (6) spineProperties:      properties of the spine change in momentum
%                           (tensor, spine elements etc.)
% (7) dragProperties:       properties of the drag force
% (8) actuTop:              vectors and matrices related to the actuation
%                           muscle at the top
% (9) actuBottom:           vectors and matrices related to the actuation
%                           muscle at the bottom
% (10) etaSol:              solution for eta to consider
% (11) etadSol:             solution for dot{eta} to consider
% (12) etadSol:             solution for ddot{eta} to consider
% (13) h:                   time step for time integration
% (14) tmax:                simulation for [0,tmax]
%
% OUTPUTS:
% (1) TI_sens:              struct containing the solutions and related
%                           information
%     
%
% Last modified: 18/05/2024, Mathieu Dubied, ETH Zurich

function TI_sens = solve_sensitivities_actu(V,p_k,PROM_Assembly,actuTop,actuBottom,etaSol,etadSol,etaddSol,h,tmax)        

    % SIMULATION PARAMETERS AND ICs _______________________________________
    s0 = zeros(size(V,2),size(p_k,1));
    sd0 = zeros(size(V,2),size(p_k,1));
    sdd0 = zeros(size(V,2),size(p_k,1));
    
    % EVALUATE PARTIAL DERIVATIVES ALONG NOMINAL SOLUTION _________________
    % actuation
%     k=400;
%     actuSignalT = @(t) k/2*(-0.2*sin(t*2*pi));    % to change below as well if needed
%     actuSignalB = @(t) k/2*(0.2*sin(t*2*pi));
%     pd_actuTop = @(a)PROM_actu_derivatives(actuTop,a);
%     pd_actuBottom = @(a)PROM_actu_derivatives(actuBottom,a);
    B1T = actuTop.B1;
    B1B = actuBottom.B1;
    B2T = actuTop.B2;
    B2B = actuBottom.B2;
    k = 1.0;
    pd_actu = @(t,q) all_derivatives_actuation_force_1(k,t,q,B1T,B2T,B1B,B2B,p_k);
    
    % TIME INTEGRATION ____________________________________________________
    
    % instantiate object for time integration
    TI_sens = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
    
    % residual function handle
    Residual_sens = @(s,sd,sdd,t)residual_linear_sens_actu(s,sd,sdd,t,PROM_Assembly, ...
        etaSol,etadSol,etaddSol, ...
        pd_actu,h);

    % time integration
    TI_sens.Integrate(s0,sd0,sdd0,tmax,Residual_sens);
    
end