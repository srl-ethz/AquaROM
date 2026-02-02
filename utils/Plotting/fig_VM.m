% ------------------------------------------------------------------------ 
% fig_VM.m
%
% Description: create a figure showing the n-th vibration mode.
%
% Last modified: 02/02/2026, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function fig = fig_VM(Mesh, nodes, elements, muscleBoundaries,  esetBC, n_VMs)

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
    

    % subplot 1: muscles' placement
    ax1 = subplot('Position',pos1);
    elementPlot = elements(:,1:4);
    L = [Lx,Ly,Lz];
    O = [-Lx,-Ly/2,-Lz/2];
    plotcube(L,O,.05,[0 0 0]);
    v1 = reshape(VMn(:,n_VMs), 3, []).';
    PlotFieldonDeformedMesh_ext(nodes, elementPlot, v1, 'factor', max(nodes(:,2)),'lineWidth',0.2);
    % ------------------------------------------------------------
    % subplot2: VM1
    ax2 = subplot('Position',pos2);

    elementPlot = elements(:,1:4);
    L = [Lx,Ly,Lz];
    O = [-Lx,-Ly/2,-Lz/2];
    plotcube(L,O,.05,[0 0 0]);
    v1 = reshape(-VMn(:,n_VMs), 3, []).';
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

    

    axis([ax1 ax2],[-0.22 0 -0.04 0.04 -0.07 0.07])

end