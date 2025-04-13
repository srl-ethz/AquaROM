% ------------------------------------------------------------------------ 
% fig_muscle_placement_VM.m
%
% Description: create a figure showing the muscle placement and the first
% vibration mode (VM)
% Synthax:
% [Mesh, nodes, elements, nsetForBC] = create_mesh(filename, myElementConstructor, propRigid)
%
% Description: Creates a FE mesh in Matlab based on an Abaqus input input
% file. Describes a set of nodes used for the boundary conditions (BC)
%
% INPUTS: 
% (1) filename:             name of Abaqus input file              
% (2) myElementConstructor: element constructor type (only tested with Tet4)
% (3) propRigid:            proporion of the fish which is rigid (used for
%                           the boundary conditions)
%
% OUTPUTS:   
% (1) mesh:             mesh object
% (2) nodes:            matrix containing the positions of the nodes
% (3) elements:         matrix containing the nodes' ID of each element
% (4) muscleBoudnaries: boundaries of the muscles, along the x axis
% (5) esetBC:           binary vector of size n_elements, 1 (0): element (not) 
%                       subject to boundary conditions (element set)
%
% Additional notes: this function is not written in the most generic way.
% It purpose is to facilitate the creation of a mesh in our examples.
%
% Last modified: 16/03/2025, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function fig = fig_muscle_placement_VM(Mesh, nodes, elements, muscleBoundaries,  esetBC)

    [Lx, Ly, Lz] = mesh_dimensions_3D(nodes);

    % MUSCLES _____________________________________________________________
    % left muscle (y>0)
    nel = size(elements,1);
    leftMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if elementCenterY>0.00 &&  elementCenterX<-Lx*muscleBoundaries(2) && elementCenterX>-Lx*muscleBoundaries(1)
            leftMuscle(el) = 1;
        end    
    end

    % right muscle (y<0)
    rightMuscle = zeros(nel,1);
    for el=1:nel
        elementCenterY = (nodes(elements(el,1),2)+nodes(elements(el,2),2)+nodes(elements(el,3),2)+nodes(elements(el,4),2))/4;
        elementCenterX = (nodes(elements(el,1),1)+nodes(elements(el,2),1)+nodes(elements(el,3),1)+nodes(elements(el,4),1))/4;
        if elementCenterY<0.00 &&  elementCenterX<-Lx*muscleBoundaries(2) && elementCenterX>-Lx*muscleBoundaries(1)
            rightMuscle(el) = 1;
        end    
    end

    % VIBRATION MODES _________________________________________________________
    % get first vibration mode
    NominalAssemblyForPlot = Assembly(Mesh);
    Mn = NominalAssemblyForPlot.mass_matrix();
    nNodes = size(nodes,1);
    u0 = zeros(Mesh.nDOFs, 1);
    [Kn,~] = NominalAssemblyForPlot.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssemblyForPlot.DATA.K = Kn;
    NominalAssemblyForPlot.DATA.M = Mn;

    % vibration modes
    n_VMs = 1;
    Kc = NominalAssemblyForPlot.constrain_matrix(Kn);
    Mc = NominalAssemblyForPlot.constrain_matrix(Mn);
    [VMn,om] = eigs(Kc, Mc, n_VMs, 'SM');
    [~,ind] = sort(sqrt(diag(om))/2/pi);
    VMn = VMn(:,ind);
    for ii = 1:n_VMs
        VMn(:,ii) = VMn(:,ii)/max(sqrt(sum(VMn(:,ii).^2,2)));
    end
    VMn = NominalAssemblyForPlot.unconstrain_vector(VMn);

    % FIGURE __________________________________________________________________
    fig = figure('units','centimeters','position',[3 3 9 5.0]);
    pos1 = [0.0,0,0.43,1];
    pos2 = [0.43,0,0.43,1];
    fixedElements = esetBC;
    
    % axis in first subplot (before the plot itself so that the arrows are
    % below the structure)
    % Create textbox
    annotation(fig,'textbox',...
        [0.445717592592594 0.621540123456792 0.0823148148148148 0.15679012345679],...
        'String',{'x'},...
        'Interpreter','latex',...
        'FitBoxToText','off',...
        'EdgeColor','none'); 
    % Create textbox
    annotation(fig,'textbox',...
        [0.251689814814815 0.801260802469136 0.0823148148148148 0.15679012345679],...
        'String',{'y'},...
        'Interpreter','latex',...
        'FitBoxToText','off',...
        'EdgeColor','none');
    % Create textbox
    annotation(fig,'textbox',...
        [0.384961419753087 0.828111111111112 0.0823148148148147 0.15679012345679],...
        'String',{'z'},...
        'Interpreter','latex',...
        'FitBoxToText','off',...
        'EdgeColor','none');

    % Create arrow
    annotation(fig,'arrow',[0.40945987654321 0.478055555555556],...
        [0.597680555555556 0.652361111111111],  'HeadSize', 8);
    % Create arrow
    annotation(fig,'arrow',[0.295787037037037 0.236010802469136],...
        [0.768777777777778 0.858736111111111],  'HeadSize', 8);
    % Create arrow
    annotation(fig,'arrow',[0.376141975308642 0.375162037037037],...
        [0.782888888888889 0.937296296296296], 'HeadSize', 8);

    % subplot 1: muscles' placement
    ax1 = subplot('Position',pos1);

    PlotMeshWith3Sets(nodes,elements, ...
        leftMuscle,'green',rightMuscle,'blue', ...
        fixedElements,'red');
    % ------------------------------------------------------------
    % subplot2: VM1
    ax2 = subplot('Position',pos2);

    elementPlot = elements(:,1:4);
    L = [Lx,Ly,Lz];
    O = [-Lx,-Ly/2,-Lz/2];
    plotcube(L,O,.05,[0 0 0]);
    v1 = reshape(-VMn(:,1), 3, []).';
    PlotFieldonDeformedMesh_ext(nodes, elementPlot, v1, 'factor', max(nodes(:,2)),'lineWidth',0.2);

    
    % Colorbar
    text('String', 'deformation', ...
    'Position', [0.085, 0,-0.118], ...  
    'Rotation', 90, ...
    'Interpreter', 'latex', ...
    'Visible','on');
    view(ax2,[-37.5 30]);
    hold(ax2,'off');
    % Set the remaining axes properties
    set(ax2,'DataAspectRatio',[1 1 1]);
    % Create colorbar
    c = colorbar(ax2,'Position',[0.888645445641528 0.25 0.0211361131025105 0.6],...
        'Limits',[0 1.05218383906483]);
    c.Ticks = [];
    c.TickLabels = {};
    % Create textbox
    annotation(fig,'textbox',...
        [0.869756121449554 0.873042527954369 0.0577659157688553 0.125944584382873],...
        'String','$+$',...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    % Create textbox
    annotation(fig,'textbox',...
        [0.869756121449556 0.154398731224928 0.0577659157688553 0.125944584382872],...
        'String','$-$',...
        'Interpreter','latex',...
        'HorizontalAlignment','center',...
        'FontSize',15,...
        'FitBoxToText','off',...
        'EdgeColor','none');

    
    % Create textarrow
    annotation(fig,'textarrow',[0.350663580246914 0.312619149103396],...
    [0.254475308641975 0.40901650496004],'String',{'rigid'},'LineWidth',1,...
    'Interpreter','latex',...
    'HorizontalAlignment','left',...
    'HeadStyle','none');
    % Create textarrow
    annotation(fig,'textarrow',[0.0958796296296296 0.136458333333333],...
    [0.800527777777778 0.634138888888889],'String',{'muscle 1'},'LineWidth',1,...
    'Interpreter','latex',...
    'HeadStyle','none');
    % Create textarrow
    annotation(fig,'textarrow',[0.198773148148148 0.198996913580247],...
    [0.129043209876543 0.331347222222222],'String',{'muscle 2'},'LineWidth',1,...
    'Interpreter','latex',...
    'HorizontalAlignment','left',...
    'HeadStyle','none');

    % Dimension: x
    text('String', '20 cm', ...
    'Position', [-0.15, 0, 0.08], ...  
    'Rotation', 20, ...
    'Interpreter', 'latex', ...
    'Visible','on');
    % Dimension: y
    text('String', '4 cm', ...
    'Position', [-0.015, 0, 0.08], ...  
    'Rotation', -35, ...
    'Interpreter', 'latex', ...
    'Visible','on');

%     % Dimension: z
%     text('String', '10 cm', ...
%     'Position', [0.03, 0.00, 0.01], ...  
%     'Rotation', -90, ...
%     'Interpreter', 'latex', ...
%     'Visible','on');
    
    % Dimension: z
    text('String', '10 cm', ...
    'Position', [-0.23, 0.00, -0.01], ...  
    'Rotation', 90, ...
    'Interpreter', 'latex', ...
    'Visible','on');

    axis([ax1 ax2],[-0.22 0 -0.04 0.04 -0.055 0.055])

end