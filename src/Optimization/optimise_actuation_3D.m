% optimise_actuation_3D
%
% Synthax:
% [xiStar,xiEvo,LrEvo] = optimise_actuation_3D(myElementConstructor,nset,nodes,elements,d,h,tmax,A,b,varargin)
%
% Description: Implementation of the optimization pipeline presented in
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
% (4) LwoBEvo:      evolution of the cost function values without barriers
%     
%
% Last modified: 03/12/2023, Mathieu Dubied, ETH Zurich

function [xiStar,pEvo,LEvo,LwoBEvo] = optimise_actuation_3D(myElementConstructor,nset,nodes,elements,d,h,tmax,A,b,varargin)

    % parse input
    [maxIteration,convCrit,convCritCost,barrierParam,gradientWeights, ...
        gStepSize,nResolve,resolveThreshold,~,~,USEJULIA] = parse_inputs(varargin{:});
    
    % NOMINAL SOLUTION ____________________________________________________
    fprintf('**************************************\n')
    fprintf('Solving for the nominal actuation...\n')
    fprintf('**************************************\n')

    %p_k = ones(size(A,2),1);
    % p_k = [1;0;0];    % actuation_force_2
    % p_k = [0.8;2*pi;0];   % actuation_force_3
%      p_k = [1;0;0;1];    % actuation_force_4
        
    p_k = [1];    % actuation_force_6, AO1

%     % actuation_force_5:
%     nParam = size(A,2);
%     p_k = 0.5*ones(nParam,1);
%     for i = 1:nParam
%         p_k(i,1) = sin((i-1)*h*2*pi);
%     end


    pResolve_k = p_k;
    deltaP_k = p_k - pResolve_k;
    pEvo = p_k;
    nParam = length(p_k);
    
    % optimisation characteristics of each constraint/parameter
    if length(barrierParam) ~= size(b,1)
        barrierParam = barrierParam*ones(1,size(b,1));
    end
    
    if length(gradientWeights) ~= nParam
        gradientWeights = gradientWeights*ones(1,nParam);
    end
   
    % Mesh
            
    MeshNominal = Mesh(nodes);
    MeshNominal.create_elements_table(elements,myElementConstructor);
    
    Lx = abs(max(nodes(:,1))-min(nodes(:,1)));  % horizontal length of airfoil
    Ly = abs(max(nodes(:,2))-min(nodes(:,2)));  % vertical length of airfoil
    Lz = abs(max(nodes(:,3))-min(nodes(:,3)));  % vertical length of airfoil
    [y_thinFish,z_smallFish,z_tail,z_head,z_linLongTail, z_notch,...
    y_tail,y_head,y_linLongTail,y_ellipseFish, x_concaveTail] = ...
    shape_variations_3D(nodes,Lx,Ly,Lz);
    U = [z_tail,z_head,y_thinFish];
 
    for l=1:length(nset)
        MeshNominal.set_essential_boundary_condition([nset{l}],1:3,0)   
    end

    % build PROM
    fprintf('____________________\n')
    fprintf('Building ROM ... \n')

    [V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom] = ...
        build_PROM_3D(MeshNominal,nodes,elements,U,USEJULIA,1,'N1t'); 
%     build_ROM_3D(MeshNominal,nodes,elements,USEJULIA);    
    

    % Solve EoMs
%     tic 
%     fprintf('____________________\n')
%     fprintf('Solving EoMs and sensitivities ...\n')
%     TI_NL_PROM = solve_EoMs_and_sensitivities_actu(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,p_k);                        
%     toc
    tic 
    fprintf('____________________\n')
    fprintf('Solving EoMs ...\n')
    TI_NL_PROM = solve_EoMs_custom_actu(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,p_k); 
    toc
    


    % Displaying solution
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
    etad = TI_NL_PROM.Solution.qd;
    etadd = TI_NL_PROM.Solution.qdd;
    
    % Solve sensitivities
    tic 
    fprintf('____________________\n')
    fprintf('Solving Sensitivities ...\n')
    TI_NL_PROM_sens = solve_sensitivities_actu(V,p_k,PROM_Assembly,actuTop,actuBottom,eta,etad,etadd,h,tmax);
    toc

%     S = TI_NL_PROM.Solution.s;
%     Sd = TI_NL_PROM.Solution.sd;
%     Sdd = TI_NL_PROM.Solution.sdd;
    S = TI_NL_PROM_sens.Solution.q;

    
    eta_0k = TI_NL_PROM.Solution.q;
    eta_k = eta;
    etad_k = etad;
    etadd_k = etadd;

    % nodes = PROM_Assembly.Mesh.nodes;
    x0 = reshape(nodes.',[],1);
    
    % computing initial cost function value
    fprintf('____________________\n')
    fprintf('Computing cost function...\n')             
    N = size(eta,2);
    dr = reduced_constant_vector(d,V,3);                
    [L,LwoB] = reduced_cost_function_w_constraints_TET4(N,eta_k,p_k,A,b,barrierParam,V);  
    LEvo = L;
    LwoBEvo = LwoB;
    nablaEvo = zeros(size(A,2),1);
    lastResolve = 0;

    for k = 1:maxIteration
        fprintf('**************************************\n')
        fprintf('Optimization loop iteration: k= %d\n',k)
        fprintf('**************************************\n')

        % possible resolve
        if check_cond_resolve(k,lastResolve,nResolve,deltaP_k,resolveThreshold,maxIteration)
            lastResolve = k;      
            pResolve_k = p_k;
                                                            
            % solve EoMs to get updated nominal solutions eta and dot{eta} (on the deformed mesh
            tic 
            fprintf('____________________\n')
            fprintf('Solving EoMs and sensitivity...\n') 
            TI_NL_PROM = solve_EoMs_and_sensitivities_actu(V,PROM_Assembly,tensors_PROM,tailProperties,spineProperties,dragProperties,actuTop,actuBottom,h,tmax,p_k);                        
            toc
            
            % get solutions
            eta_0k = TI_NL_PROM.Solution.q;
            eta_k = TI_NL_PROM.Solution.q;
            etad_k = TI_NL_PROM.Solution.qd;
            etadd_k = TI_NL_PROM.Solution.qdd;
            S = TI_NL_PROM.Solution.s;
            Sd = TI_NL_PROM.Solution.sd;
            Sdd = TI_NL_PROM.Solution.sdd;
    
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
                % plot solution
                uTail = zeros(3,tmax/h);
                for a=1:tmax/h
                    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*eta_k(:,a);
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
                eta_k = eta_0k + S*deltaP_k;
                % plot solution
                uTail = zeros(3,tmax/h);
                for a=1:tmax/h
                    uTail(:,a) = V(tailProperties.tailNode*3-2:tailProperties.tailNode*3,:)*eta_k(:,a);
                end 
                subplot(2,1,1);
                plot(timePlot,x0Tail+uTail(1,:),'DisplayName',strcat('k=',num2str(k)))
                legend
                drawnow

                subplot(2,1,2);
                plot(timePlot,uTail(2,:),'DisplayName',strcat('k=',num2str(k)))
                legend
                drawnow
            end
            
        end 

        % compute cost function and its gradient
        fprintf('____________________\n')
        fprintf('Computing cost function and its gradient...\n')   
        nablaLr = gradient_cost_function_w_constraints_TET4(p_k,eta_k,S,A,b,barrierParam,V);
        [L,LwoB] = reduced_cost_function_w_constraints_TET4(N,eta_k,p_k,A,b,barrierParam,V);
        LEvo = [LEvo, L];
        LwoBEvo = [LwoBEvo, LwoB];

        % specific for AO2
        %nablaLr(1:int32(tmax/h*0.2)) = 0;

        nablaEvo = [nablaEvo,nablaLr];

        % update optimal parameter
        fprintf('____________________\n')
        fprintf('Updating optimal parameter...\n') 
        gradientWeights = adapt_learning_rate(nablaEvo,gradientWeights);

        p_k = p_k - gStepSize*diag(gradientWeights)*nablaLr

        %p_k = clip_infeasible_parameters(p_k,A,b);
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

        deltaP_k = p_k - pResolve_k;    
        pEvo = [pEvo,p_k];
        
        % possible exit conditions
        if size(p_k,1) >1
            if norm(pEvo(:,end)-pEvo(:,end-1))<convCrit
                fprintf('Convergence criterion of %.3f (parameters) fulfilled\n',convCrit)
                break
            elseif length(LEvo)>5
                if norm(LEvo(end) - mean(LEvo(end-4:end))) < convCritCost
                    fprintf('Convergence criterion of %.3f (cost) fulfilled\n',convCritCost)
                    break
                end
            end
        else
            if norm(pEvo(end)-pEvo(end-1))<convCrit
                fprintf('Convergence criterion of %.3f (parameters) fulfilled\n',convCrit)
                break
            elseif length(LEvo)>5
                if norm(LEvo(end) - mean(LEvo(end-4:end))) < convCritCost
                    fprintf('Convergence criterion of %.3f (cost) fulfilled\n',convCritCost)
                    break
                end
            end
        end

        if k == maxIteration
            fprintf('Maximum number of %d iterations reached\n',maxIteration)
        end
    end
    
    xiStar = p_k;

end

% Parse input _____________________________________________________________
function [maxIteration,convCrit,convCritCost,barrierParam,gradientWeights,gStepSize,nResolve,resolveThreshold,FORMULATION,VOLUME,USEJULIA] = parse_inputs(varargin)
    defaultMaxIteration = 50;
    defaultConvCrit = 0.001;
    defaultConvCritCost = 0.1;
    defaultBarrierParam = 500;
    defaultGradientWeights = 1;
    defaultGStepSize = 0.1;
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
    nResolve = p.Results.nResolve;
    resolveThreshold = p.Results.resolveThreshold;
    FORMULATION = p.Results.FORMULATION;
    VOLUME = p.Results.VOLUME;
    USEJULIA = p.Results.USEJULIA;
end

% Check condition for rebuild _____________________________________________
function cond = check_cond_resolve(k,lastRebuild,nRebuild, xiRebuild_k, ...
                                    rebuildThreshold,maxIteration)
    cond = 0;

    if mod(k-lastRebuild,nRebuild) == 0 
        cond = 1;
        fprintf('____________________\n')
        fprintf('Resolving EoMs (max lin. iterations) ...\n')
    elseif any(xiRebuild_k > rebuildThreshold)
        cond = 1;
        fprintf('____________________\n')
        fprintf('Resolving EoMs (xi>threshold) ...\n')
    elseif maxIteration-k<0.25*maxIteration ...
            && mod(k-lastRebuild,int16(nRebuild/2)) == 0
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
            fprintf('Adapting learning rate (change in gradient sign) - xi%d...\n',p)
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





