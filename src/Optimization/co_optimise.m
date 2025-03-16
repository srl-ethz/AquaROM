
% co_optimise
%
% Synthax:
% [xiStar,xiEvo,LEvo,LwoBEvo] = co_optimise(myElementConstructor,nset,nodes,elements,U,h,tmax,A,b,varargin)
%
% Description: Implementation of the (shape and actuation) co-optimisation 
% pipeline presented in the paper
%
% INPUTS: 
% (1) myElementConstructor: defines the type of element and material
%                           properties
% (2) nset:                 set of element to apply boundary conditions             
% (3) nodes:                nodes and their coordinates
% (4) elements:             elements described by their nodes
% (5) U:                    shape variation basis
% (6) h:                    time step for time integration
% (7) tmax:                 simulation for [0,tmax]
% (8)-(9) A,b               constraints on xi of the form Axi<b
% (10) nPShape              number of shape parameters
% (11) nPActu               number of actuation parameters
%
% possible additional name-value pair arguments
% (12) maxIteration:maximum number of iterations
% (13) convCrit:    convergence criterium. Norm between two successive
%                   optimal paramter vectors
% (14) FORMULATION: order of the Neumann approximation (N0/N1/N1t)
% (15) VOLUME:      integration over defected (1) or nominal volume (0)
% (16) USEJULIA:    use of JULIA (1) for the computation of internal forces
%                   tensors
% (17) barrierParam:parameter to scale the barrier function for the 
%                   constraints (1/barrierParam)
% (18) gStepSize:   step size used in the gradient descent algorithm
% (19) nRebuild:    number of step between each re-build of a PROM
%
% OUTPUTS:
% (1) xiStar:       optimal shape parameter(s) (scalar or vector)
% (2) xiEvo:        evolution of the optimal shape parameter(s)
% (3) LrEvo:        evolution of the cost function values
%     
%
% Last modified: 07/01/2024, Mathieu Dubied, ETH Zurich

function [pStar,pEvo,LEvo,LwoBEvo] = co_optimise(myElementConstructor,nset,nodes,elements,U,h,tmax,A,b,nPShape,nPActu,varargin)

    % parse input
    [maxIteration,convCrit,convCritCost,barrierParam,gradientWeights, ...
        gStepSize,nRebuild,...
        rebuildThreshold,nResolve,resolveThreshold,...
        FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin{:});
    rebuildThresholdSwitch = 0;
    
    % NOMINAL SOLUTION ____________________________________________________
    fprintf('**************************************\n')
    fprintf('Solving for the nominal structure...\n')
    fprintf('**************************************\n')

    % shape parameters
    xi_k = zeros(size(U,2),1);%[-0.0694;-0.487;0.2420;-0.287;0.058;0.3603;0.475;0.2852];%zeros(size(U,2),1);
%     xi_k(8) = 0.2;%[-0.3941;0.3941]
    xiAtLastRebuild = xi_k;%zeros(size(U,2),1);
    deltaXi_k = xi_k - xiAtLastRebuild;
    
    % actuation parameters
    pActu_k = [0.2;2;0.5];%[0.2517;2.3089;0.6876];%[0.2;2;0.1];% [1.0240;0.0852;0.1689;1.2457];%[1.1;0.0;0.18;1.20];%[1;0;0;1];    % actuation_force_4
    pActuAtLastResolve = pActu_k;
    deltaPActu_k = pActu_k - pActuAtLastResolve;
    
    % merging parameters
    nParam = nPShape + nPActu;
    p_k = [xi_k; pActu_k];
    deltaP_k = [deltaXi_k;deltaPActu_k];
    pEvo = p_k;
    
    % optimisation characteristics of each constraint/parameter
    if length(barrierParam) ~= size(b)
        barrierParam = barrierParam*ones(1,nParam);
    end
    
    if length(gradientWeights) ~= nParam
        gradientWeights = gradientWeights*ones(1,nParam);
    end
    
    % Mesh      
    MeshNominal = Mesh(nodes);
    MeshNominal.create_elements_table(elements,myElementConstructor);
 
    for l=1:length(nset)
        MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)   
    end
    
    % to test
    df = U*xi_k;                       % displacement fields introduced by defects
    ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
    nodes_defected = nodes + ddf;    % nominal + d ---> defected 
    svMesh = Mesh(nodes_defected);
    svMesh.create_elements_table(elements,myElementConstructor);
    for l=1:length(nset)
        svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
    end

    % build PROM
    fprintf('____________________\n')
    fprintf('Building PROM ... \n')

    [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
    build_PROM_3D(svMesh,nodes_defected,elements,U,USEJULIA,VOLUME,FORMULATION);      
    

    % Solve EoMs
    tic 
    fprintf('____________________\n')
    fprintf('Solving EoMs and sensitivities...\n') 
    TI_NL_PROM = solve_EoMs_and_sensitivities_co(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,pActu_k); 
                
    toc

    uTail = zeros(3,tmax/h);
    timePlot = linspace(0,tmax-h,tmax/h);
    x0Tail = min(nodes(:,1));
    for a=1:tmax/h
        uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
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


    
    % Retrieving solutions    
    eta = TI_NL_PROM.Solution.q;
    S = TI_NL_PROM.Solution.s;
    eta_0k = TI_NL_PROM.Solution.q;
    eta_k = eta;
    
    % computing initial cost function value
    fprintf('____________________\n')
    fprintf('Computing cost function...\n') 
    N = size(eta,2);
    [L,LwoB] = reduced_cost_function_w_constraints_TET4(N,eta_k,p_k,A,b,barrierParam,V);  
    LEvo = L;
    LwoBEvo = LwoB;
    nablaEvo = zeros(size(A,2),1);

    lastRebuild = 0;
    lastResolve = 0;
    maxItACTIVE = 1;

    for k = 1:maxIteration
        fprintf('**************************************\n')
        fprintf('Optimization loop iteration: k= %d\n',k)
        fprintf('**************************************\n')

        % possible rebuilding of a PROM
        if check_cond_rebuild(k,lastRebuild,nRebuild,deltaXi_k,rebuildThreshold,maxIteration)
            lastRebuild = k;
            lastResolve = k;
            maxItACTIVE = 0;
         
            % update defected mesh nodes
            df = U*xi_k;                       % displacement fields introduced by defects
            ddf = [df(1:3:end) df(2:3:end) df(3:3:end)]; 
            nodes_defected = nodes + ddf;    % nominal + d ---> defected 
            svMesh = Mesh(nodes_defected);
            svMesh.create_elements_table(elements,myElementConstructor);
            for l=1:length(nset)
                svMesh.set_essential_boundary_condition([nset{l}],1:3,0)   
            end

            % build PROM
            [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
                 build_PROM_3D(svMesh,nodes_defected,elements,U,USEJULIA,VOLUME,FORMULATION);
            
            xiAtLastRebuild = xi_k;      
              
            % solve EoMs to get updated nominal solutions eta and dot{eta} (on the deformed mesh
            tic 
            fprintf('____________________\n')
            fprintf('Solving EoMs and sensitivity...\n') 
            TI_NL_PROM = solve_EoMs_and_sensitivities_co(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,pActu_k);                        
            toc
             
            eta_0k = TI_NL_PROM.Solution.q;
            eta_k = TI_NL_PROM.Solution.q;
            S = TI_NL_PROM.Solution.s;
            
            pActuAtLastResolve = p_k(nPShape+1:end);
    
            N = size(eta_k,2);

            uTail = zeros(3,tmax/h);
            for a=1:tmax/h
                uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
            end 
            subplot(2,1,1);
            plot(timePlot,x0Tail+uTail(1,:),'DisplayName',strcat('k=',num2str(k)))
            legend
            drawnow
     
            subplot(2,1,2);
            plot(timePlot,uTail(2,:),'DisplayName',strcat('k=',num2str(k)))
            legend
            drawnow
          
        % possible resolve of the EoMs (without rebuilding the PROM)
        elseif check_cond_resolve(k,lastResolve,nResolve,deltaPActu_k,resolveThreshold,maxIteration,maxItACTIVE)
            lastResolve = k;      
            pActuAtLastResolve = p_k(nPShape+1:end);
                                                            
            % solve EoMs to get updated nominal solutions eta and dot{eta} (on the deformed mesh
            tic 
            fprintf('____________________\n')
            fprintf('Solving EoMs and sensitivity...\n') 
            TI_NL_PROM = solve_EoMs_and_sensitivities_co(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,pActu_k);                        
            toc
            
            % get solutions
            eta_0k = TI_NL_PROM.Solution.q;
            eta_k = TI_NL_PROM.Solution.q;
            S = TI_NL_PROM.Solution.s;
    
            N = size(eta_k,2);

            % plot solution
            uTail = zeros(3,tmax/h);
            for a=1:tmax/h
                uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*TI_NL_PROM.Solution.q(:,a);
            end 
            subplot(2,1,1);
            plot(timePlot,x0Tail+uTail(1,:),'DisplayName',strcat('k=',num2str(k)))
            legend
            drawnow
     
            subplot(2,1,2);
            plot(timePlot,uTail(2,:),'DisplayName',strcat('k=',num2str(k)))
            legend
            drawnow
        else
            % approximate new solution under new xi, using sensitivity
            fprintf('____________________\n')
            fprintf('Approximating solutions...\n')
            if size(p_k,1)>1
                S=tensor(S);
                eta_k = eta_0k + double(ttv(S,deltaP_k,2));
%                 uTail = zeros(3,tmax/h);
%                 for a=1:tmax/h
%                     uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*eta_k(:,a);
%                 end 
%                 subplot(2,1,1);
%                 plot(timePlot,x0Tail+uTail(1,:),'DisplayName',strcat('k=',num2str(k)))
%                 legend
%                 drawnow
% 
%                 subplot(2,1,2);
%                 plot(timePlot,uTail(2,:),'DisplayName',strcat('k=',num2str(k)))
%                 legend
%                 drawnow
            else
                eta_k = eta_0k + S*deltaP_k;
            end
            
        end 

        % compute cost function and its gradient
        fprintf('____________________\n')
        fprintf('Computing cost function and its gradient...\n')   
        nablaLr = gradient_cost_function_w_constraints_TET4(p_k,eta_k,S,A,b,barrierParam,V);
        [L,LwoB] = reduced_cost_function_w_constraints_TET4(N,eta_k,p_k,A,b,barrierParam,V);

        LEvo = [LEvo, L];
        LwoBEvo = [LwoBEvo, LwoB];
        nablaEvo = [nablaEvo,nablaLr];

        % update optimal parameters
        fprintf('____________________\n')
        fprintf('Updating optimal parameter...\n') 
        updatedGradientWeights = adapt_learning_rate(nablaEvo,gradientWeights);

        if ~all(gradientWeights == updatedGradientWeights)  
%             if  rebuildThresholdSwitch ==0 && k>20
%                 rebuildThreshold = rebuildThreshold/2;
%                 rebuildThresholdSwitch = 1;
%             end
            gradientWeights = updatedGradientWeights;
            
        end
        
        % update parameter
        updateVector = gStepSize*diag(gradientWeights)*nablaLr;
        p_k = p_k - updateVector
               
        p_k_clipped = clip_infeasible_parameters(p_k,A,b);
        if ~all(p_k_clipped == p_k)
%             xiAtLastRebuild = xiAtLastRebuild + (p_k_clipped(1:nPShape) - xi_k);
%             if ~all(p_k_clipped(1:nPShape) == p_k(1:nPShape))
%                 %rebuildThreshold = rebuildThreshold/2;
%             end
           idxToChange = find(p_k~=p_k_clipped);
           for idx = 1:length(idxToChange)
               gradientWeights(idxToChange(idx)) = ...
                   0.5*gradientWeights(idxToChange(idx));
               fprintf('Adapting learning rate for xi%d to %.3f...\n',...
                   idxToChange(idx),gradientWeights(idxToChange(idx)))
               if nParam*2 == length(b)
                    barrierParam(idxToChange(idx)*2-1:idxToChange(idx)*2) = ...
                        0.5*barrierParam(idxToChange(idx)*2-1:idxToChange(idx)*2);
                    fprintf('Decreasing barrier parameter %d to %.3f...\n',...
                        idxToChange(idx),barrierParam(idxToChange(idx)*2-1))
               end
               
           end
            
            p_k = p_k_clipped;
        end
         
        % update shape parameter
        xi_k = p_k(1:nPShape);
        deltaXi_k = xi_k - xiAtLastRebuild;
        
        % update actuation parameters
        pActu_k = p_k(nPShape+1:end);
        deltaPActu_k = pActu_k - pActuAtLastResolve;
        
        % difference with the last rebuild/resolve
        deltaP_k = [deltaXi_k;deltaPActu_k];
        
        % overall evolution
        pEvo = [pEvo,p_k];
        maxItACTIVE = 1;
        
        % possible exit conditions
        if size(p_k,1) >1
            if norm(pEvo(:,end)-pEvo(:,end-1))<convCrit
                fprintf('Convergence criterion of %.3f (parameters) fulfilled\n',convCrit)
                break
            elseif length(LEvo)>7
                if var(LEvo(end-6:end)) < convCritCost
                    fprintf('Convergence criterion of %.3f (cost) fulfilled\n',convCritCost)
                    break
                end
            end
        else
            if norm(pEvo(end)-pEvo(end-1))<convCrit
                fprintf('Convergence criterion of %.3f(parameters) fulfilled\n',convCrit)
                break
            elseif length(LEvo)>7
                if var(LEvo(end-6:end)) < convCritCost
                    fprintf('Convergence criterion of %.3f (cost) fulfilled\n',convCritCost)
                    break
                end
            end
        end

        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end
    
    pStar = p_k;

end

% Parse input _____________________________________________________________
function [maxIteration,convCrit,convCritCost,barrierParam, ...
    gradientWeights,gStepSize,nRebuild,rebuildThreshold, nResolve, ...
    resolveThreshold, FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin)
    defaultMaxIteration = 50;
    defaultConvCrit = 0.001;
    defaultConvCritCost = 0.1;
    defaultBarrierParam = 500;
    defaultGradientWeights = 1;
    defaultGStepSize = 0.1;
    defaultNRebuild = 10;
    defaultRebuildThreshold = 0.2;
    defaultNResolve = 10;
    defaultResolveThreshold = 0.2;
    defaultFORMULATION = 'N1';
    defaultVOLUME = 1;
    defaultUSEJULIA = 0; 
    p = inputParser;
    addParameter(p,'maxIteration',defaultMaxIteration, @(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','integer','positive'}) );
    addParameter(p,'convCrit',defaultConvCrit,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'convCritCost',defaultConvCritCost,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'barrierParam',defaultBarrierParam,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'gradientWeights',defaultGradientWeights,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'gStepSize',defaultGStepSize,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'nRebuild',defaultNRebuild,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'rebuildThreshold',defaultRebuildThreshold,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'nResolve',defaultNResolve,@(x)validateattributes(x, ...
                    {'numeric'},{'nonempty','positive'}) );
    addParameter(p,'resolveThreshold',defaultResolveThreshold,@(x)validateattributes(x, ...
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
    convCritCost = p.Results.convCritCost;
    barrierParam = p.Results.barrierParam;
    gradientWeights = p.Results.gradientWeights;
    gStepSize = p.Results.gStepSize;
    nRebuild = p.Results.nRebuild;
    rebuildThreshold = p.Results.rebuildThreshold;
    nResolve = p.Results.nResolve;
    resolveThreshold = p.Results.resolveThreshold;
    FORMULATION = p.Results.FORMULATION;
    VOLUME = p.Results.VOLUME;
    USEJULIA = p.Results.USEJULIA;
end

% Check conditions for rebuild ____________________________________________
function cond = check_cond_rebuild(k,lastRebuild,nRebuild, deltaXi_k, ...
                                    rebuildThreshold,maxIteration)
    cond = 0;

    if mod(k-lastRebuild,nRebuild) == 0 
        cond = 1;
        fprintf('____________________\n')
        fprintf('Rebuilding PROM (max lin. iterations) ...\n')
    elseif any(abs(deltaXi_k) > rebuildThreshold)
        cond = 1;
        criticalParams = find(abs(deltaXi_k) > rebuildThreshold);
        fprintf('____________________\n')
        fprintf('Rebuilding PROM (delta xi>threshold) for \n') 
        fprintf(' xi%.0f \n', criticalParams) 
    elseif maxIteration-k<0.2*maxIteration ...
            && mod(k-lastRebuild,int16(nRebuild/1.33)) == 0
        cond = 1;
        fprintf('____________________\n')
        fprintf('Rebuilding PROM (max lin. iteration - close to end) ...\n')
    end
    
end

% Check conditions for resolve ____________________________________________
function cond = check_cond_resolve(k,lastResolve,nResolve, deltaPActu_k, ...
                                    resolveThreshold,maxIteration, ...
                                    maxItACTIVE)
    cond = 0;

    if mod(k-lastResolve,nResolve) == 0 
        cond = 1;
        fprintf('____________________\n')
        fprintf('Resolving EoMs (max lin. iterations) ...\n')
    elseif any(abs(deltaPActu_k) > resolveThreshold)
        cond = 1;
        fprintf('____________________\n')
        fprintf('Resolving EoMs(delta pActu>threshold) ...\n')
    elseif maxIteration-k<0.2*maxIteration ...
            && mod(k-lastResolve,int16(nResolve/1.33)) == 0 ...
            && maxItACTIVE == 1
        cond = 1;
        fprintf('____________________\n')
        fprintf('Resolving EoMs (max lin. iteration - close to end) ...\n')
    end
    
end


% Adapt gradient __________________________________________________________
function gradientWeights = adapt_learning_rate(nablaEvo,currentGradientWeights)
    nParam = size(nablaEvo,1);
    gradientWeights = currentGradientWeights;
    for p = 1:nParam
        if sign(nablaEvo(p,end-1)) ~= sign(nablaEvo(p,end)) ...
                && nablaEvo(p,end-1) ~= 0
            fprintf('')
            gradientWeights(p) = 0.5*currentGradientWeights(p);
            fprintf('Adapting learning rate for xi%d to %.3f...\n',p,gradientWeights(p))
        end
    end
end

% Clip infeasible parameters ______________________________________________
% Note: only work for constraint containing a single parameter
function param = clip_infeasible_parameters(p,A,b)

    param = p;

    % check if problem is infeasible
    if ~all(A*p<b)
        constrIdxToClip = find(A*p>b);   
        
        % iterate over violated constraints
        for i = 1:length(constrIdxToClip)
            constrIdx = constrIdxToClip(i);

            % only consider constraints containing a single parameter
            if length(find(A(constrIdx,:))) == 1
                paramIdxToClip = find(A(constrIdx,:));
                param(paramIdxToClip) = sign(A(constrIdx,paramIdxToClip))*b(constrIdx) - sign(A(constrIdx,paramIdxToClip))*0.05*b(constrIdx);
                fprintf('Clipping  parameter %d to the value %d \n',paramIdxToClip,param(paramIdxToClip))
            end
        end
    end
 end

