% optimise_shape_2D
%
% Synthax:
% [xiStar,xiEvo,LrEvo] = optimise_shape_2D(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline 4 presented in
% the paper, based on Newton's method. 
%
% INPUTS: 
% (1) myElementConstructor: defines the type of element and material
%                           properties
% (2) nset:                 set of element to apply boundary conditions             
% (3) nodes:                nodes and their coordinates
% (4) elements:             elements described by their nodes
% (5) U:                    shape variation basis
% (6) d:                    forward swimming direction
% (7) h:                    time step for time integration
% (8) tmax:                 simulation for [0,tmax]
% (9)-(10) A,b              constraints on xi of the form Axi<b
%
% possible additional name-value pair arguments
% (11) maxIteration:maximum number of iterations
% (12) convCrit:    convergence criterium. Norm between two successive
%                   optimal paramter vectors
% (13) FORMULATION: order of the Neumann approximation (N0/N1/N1t)
% (14) VOLUME:      integration over defected (1) or nominal volume (0)
% (15) USEJULIA:    use of JULIA (1) for the computation of internal forces
%                   tensors
% (16) barrierParam:parameter to scale the barrier function for the 
%                   constraints (1/barrierParam)
% (17) gStepSize:   step size used in the gradient descent algorithm
% (18) nRebuild:    number of step between each re-build of a PROM
%
% OUTPUTS:
% (1) xiStar:       optimal shape parameter(s) (scalar or vector)
% (2) xiEvo:        evolution of the optimal shape parameter(s)
% (3) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 15/10/2023, Mathieu Dubied, ETH Zurich

function [xiStar,xiEvo,LrEvo] = optimise_shape_2D(myElementConstructor,nset,nodes,elements,U,d,h,tmax,A,b,varargin)

    % parse input
    [maxIteration,convCrit,barrierParam,gStepSize,nRebuild,FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin{:});
    
    % NOMINAL SOLUTION ____________________________________________________
    fprintf('**************************************\n')
    fprintf('Solving for the nominal structure...\n')
    fprintf('**************************************\n')

    xi_k = zeros(size(U,2),1);
    xiRebuild_k = zeros(size(U,2),1);
    xiEvo = xi_k;

    % Mesh
            
    MeshNominal = Mesh(nodes);
    MeshNominal.create_elements_table(elements,myElementConstructor);
 
    for l=1:length(nset)
        MeshNominal.set_essential_boundary_condition([nset{l}],1:2,0)   
    end

    % build PROM
    fprintf('____________________\n')
    fprintf('Building PROM ... \n')
    tic
    mTilde = 10;
    [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM(MeshNominal,nodes,elements,mTilde,U,USEJULIA,VOLUME,FORMULATION);      
    toc

    % Solve EoMs
    tic 
    fprintf('____________________\n')
    fprintf('Solving EoMs...\n') 
    TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax);                        
    toc

    uTail = zeros(2,tmax/h);
    timePlot = linspace(0,tmax-h,tmax/h);
    x0Tail = nodes(tailProperties.tailNode,1);
    for a=1:tmax/h
        uTail(:,a) = V(tailProperties.tailNode*2-1:tailProperties.tailNode*2,:)*TI_NL_PROM.Solution.q(:,a);
    end

    figure
    subplot(2,1,1);
    plot(timePlot,x0Tail+uTail(1,:),'DisplayName','k=0')
    hold on
    xlabel('Time [s]')
    ylabel('x-position tail node')
    legend('Location','northwest')

    subplot(2,1,2);
    plot(timePlot,uTail(2,:),'DisplayName','k=0')
    hold on
    xlabel('Time [s]')
    ylabel('y-position tail node')
    legend('Location','southwest')
    drawnow


    % Solve sensitivity equation 
    tic
    fprintf('____________________\n')
    fprintf('Solving sensitivity...\n') 
    TI_sens = solve_sensitivities(V,xi_k,PROM_Assembly, ...
       tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom, ...
       TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd,TI_NL_PROM.Solution.qdd, ...
       h,tmax);
    toc
    
    % Retrieving solutions    
    eta = TI_NL_PROM.Solution.q;
    etad = TI_NL_PROM.Solution.qd;
    etadd = TI_NL_PROM.Solution.qdd;
    S = TI_sens.Solution.q;
    Sd = TI_sens.Solution.qd;
    Sdd = TI_sens.Solution.qdd;
    eta_k = eta;
    etad_k = etad;
    etadd_k = etadd;

    nodes = PROM_Assembly.Mesh.nodes;
    x0 = reshape(nodes.',[],1);
    
    % computing initial cost function value
    fprintf('____________________\n')
    fprintf('Computing cost function...\n') 
    N = size(eta,2);
    dr = d;%reduced_constant_vector(d,V);      
    Lr = reduced_cost_function_w_constraints(N,tailProperties,spineProperties,x0,eta,etad,etadd,xiRebuild_k,xi_k,dr,A,b,barrierParam,V);  
    LrEvo = Lr;

    for k = 1:maxIteration
        fprintf('**************************************\n')
        fprintf('Optimization loop iteration: k= %d\n',k)
        fprintf('**************************************\n')

        % possible rebuilding of a PROM
        if mod(k,nRebuild) == 0
            fprintf('____________________\n')
            fprintf('Rebuilding PROM ... \n')

            % update defected mesh nodes
            df = U*xi_k;                       % displacement fields introduced by defects
            ddf = [df(1:2:end) df(2:2:end)]; 
            nodes_defected = nodes + ddf;    % nominal + d ---> defected 
            svMesh = Mesh(nodes_defected);
            svMesh.create_elements_table(elements,myElementConstructor);
            for l=1:length(nset)
                svMesh.set_essential_boundary_condition([nset{l}],1:2,0)   
            end

            % build PROM
            [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
                 build_PROM(svMesh,nodes_defected,elements,mTilde,U,USEJULIA,VOLUME,FORMULATION);
              
            xiRebuild_k = zeros(size(U,2),1);   % reset local xi to 0 as we rebuild the ROM
            
            dr = reduced_constant_vector(d,V);
    
            % solve EoMs to get updated nominal solutions eta and dot{eta} (on the deformed mesh
            tic 
            fprintf('____________________\n')
            fprintf('Solving EoMs...\n') 
            TI_NL_PROM = solve_EoMs(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax);
            toc
                
            eta_k = TI_NL_PROM.Solution.q;
            etad_k = TI_NL_PROM.Solution.qd;
            etadd_k = TI_NL_PROM.Solution.qdd;
    
            N = size(eta_k,2);
            
            % solve sensitivity equation 
            tic
            fprintf('____________________\n')
            fprintf('Solving sensitivity...\n') 
            TI_sens = solve_sensitivities(V,xiRebuild_k,PROM_Assembly, ...
               tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom, ...
               TI_NL_PROM.Solution.q,TI_NL_PROM.Solution.qd,TI_NL_PROM.Solution.qdd, ...
               h,tmax);
            toc 
    
            S = TI_sens.Solution.q;
            Sd = TI_sens.Solution.qd;
            Sdd = TI_sens.Solution.qdd;
        else
            % approximate new solution under new xi, using sensitivity
            tic 
            fprintf('____________________\n')
            fprintf('Approximating solutions...\n')
            if size(xi_k,1)>1
                S=tensor(S);
                eta_k = eta_k + double(ttv(S,xiRebuild_k,2));
            else
                eta_k = eta_k + S*xiRebuild_k;
            end
            toc
        end 

        % compute cost function and its gradient
        fprintf('____________________\n')
        fprintf('Computing cost function and its gradient...\n') 
        nablaLr = gradient_cost_function_w_constraints(dr,xiRebuild_k,xi_k,x0,eta_k,etad_k,etadd_k,S,Sd,Sdd,tailProperties,spineProperties,A,b,barrierParam,V);
        LrEvo = [LrEvo, reduced_cost_function_w_constraints(N,tailProperties,spineProperties,x0,eta,etad,etadd,xiRebuild_k,xi_k,dr,A,b,barrierParam,V)];
        
        % update optimal parameter
        fprintf('____________________\n')
        fprintf('Updating optimal parameter...\n') 
        xi_k = xi_k - gStepSize*nablaLr
        xiRebuild_k = xiRebuild_k - gStepSize*nablaLr;

        xiEvo = [xiEvo,xi_k];
        
        % possible exit conditions
        if size(xi_k,1) >1
            if norm(xiEvo(:,end)-xiEvo(:,end-1))<convCrit
                fprintf('Convergence criterion of %.2g fulfilled\n',convCrit)
                break
            end
        else
            if norm(xiEvo(end)-xiEvo(end-1))<convCrit
                fprintf('Convergence criterion of %.2g fulfilled\n',convCrit)
                break
            end
        end

        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end
    
    xiStar = xi_k;

end

% parse input
function [maxIteration,convCrit,barrierParam,gStepSize,nRebuild,FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin)
defaultMaxIteration = 50;
defaultConvCrit = 0.001;
defaultBarrierParam = 500;
defaultGStepSize = 0.1;
defaultNRebuild = 10;
defaultFORMULATION = 'N1';
defaultVOLUME = 1;
defaultUSEJULIA = 0; 
p = inputParser;
addParameter(p,'maxIteration',defaultMaxIteration, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'convCrit',defaultConvCrit,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'barrierParam',defaultBarrierParam,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'gStepSize',defaultGStepSize,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'nRebuild',defaultNRebuild,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
addParameter(p,'FORMULATION',defaultFORMULATION,@(x)validateattributes(x, ...
                {'char'},{'nonempty'}))
addParameter(p,'VOLUME',defaultVOLUME,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
addParameter(p,'USEJULIA',defaultUSEJULIA,@(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );

parse(p,varargin{:});

maxIteration = p.Results.maxIteration;
convCrit = p.Results.convCrit;
barrierParam = p.Results.barrierParam;
gStepSize = p.Results.gStepSize;
nRebuild = p.Results.nRebuild;
FORMULATION = p.Results.FORMULATION;
VOLUME = p.Results.VOLUME;
USEJULIA = p.Results.USEJULIA;

end