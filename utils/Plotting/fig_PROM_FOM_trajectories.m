% ------------------------------------------------------------------------ 
% fig_PROM_FOM_trajectories.m
% 
% Description: Create a figure showing the time trajectories of the tail
% oscillation (y-position) and the head hoirzontal displacement
% (x-position)
%
% INPUTS: 
% (1) nodes:        matrix containing the positions of the nodes
% (2) elements:     matrix containing the nodes' ID of each element
% (3) colorSpine:   color of the spine elements
% (4) colorTail:    color of the tail element
% (5) figsize:      position and size of the figure [x, y, width, height]
%
% OUTPUTS:   
% (1) fig: a figure
%
% Last modified: 02/05/2025, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function fig = fig_PROM_FOM_trajectories(uHeadTail_sol, Lx, Ly, figsize, y_shift_lgd)
    
    % read solutions to plot
    uHead_FOM = uHeadTail_sol.uHead_FOM;
    uHead_ROM = uHeadTail_sol.uHead_ROM;
    uHead_PROM_3 = uHeadTail_sol.uHead_PROM_3;
    uHead_PROM_5 = uHeadTail_sol.uHead_PROM_4;
    uHead_PROM_8 = uHeadTail_sol.uHead_PROM_5;
    
    uTail_FOM = uHeadTail_sol.uTail_FOM;
    uTail_ROM = uHeadTail_sol.uTail_ROM;
    uTail_PROM_3 = uHeadTail_sol.uTail_PROM_3;
    uTail_PROM_5 = uHeadTail_sol.uTail_PROM_4;
    uTail_PROM_8 = uHeadTail_sol.uTail_PROM_5;
    
    % figure with specific size
    figure = figure('units','centimeters','position',figsize);

    % x-position (Head)
    subplot(2,1,1);
    hold on;
    p1 = plot(timePlot, uHead_FOM(1,:), '--', 'DisplayName', 'FOM', 'LineWidth', 1.0);
    p2 = plot(timePlot, uHead_ROM(1,:), 'DisplayName', 'ROM', 'LineWidth', 1.0);
    p3 = plot(timePlot, uHead_PROM_3(1,:), '-.', 'DisplayName', 'PROM 3p', 'LineWidth', 1.0);
    p4 = plot(timePlot, uHead_PROM_5(1,:), '-.', 'DisplayName', 'PROM 5p', 'LineWidth', 1.0);
    p5 = plot(timePlot, uHead_PROM_8(1,:), '-.', 'DisplayName', 'PROM 8p', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    % --- Left y-axis: raw cm
    yyaxis left
    ylTop = ylabel('Head x-position [cm]');
    yt = yticks;
    yl = ylim;
    % --- Right y-axis: scaled to % of body length
    yyaxis right
    ax = gca; % get current axes
    ylim(yl);
    yticks(yt)
    yticklabels(yt /(Lx*100)*100); 
    ax.YColor = 'k';
    ylTopRight = ylabel('\% of body length', 'color', 'k');
    % --- Back to left for legend
    yyaxis left
    lgd = legend([p1,p2,p3,p4,p5],{'FOM', 'ROM',  'PROM 3p', 'PROM 5p', 'PROM 8p'},'Location','northwest', 'NumColumns',1, 'interpreter', 'latex');
    lgd.Position(2) = lgd.Position(2) + y_shift_lgd;

    % y-position (Tail)
    subplot(2,1,2);
    hold on;
    plot(timePlot, uTail_FOM(2,:), '--', 'DisplayName', 'FOM', 'LineWidth', 1.0);
    plot(timePlot, uTail_ROM(2,:), 'DisplayName', 'ROM', 'LineWidth', 1.0);
    plot(timePlot, uTail_PROM_3(2,:), '-.', 'DisplayName', 'PROM 3p', 'LineWidth', 1.0);
    plot(timePlot, uTail_PROM_5(2,:), '-.', 'DisplayName', 'PROM 5p', 'LineWidth', 1.0);
    plot(timePlot, uTail_PROM_8(2,:), '-.', 'DisplayName', 'PROM 8p', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    % --- Left y-axis: raw cm
    yyaxis left
    ylBottom = ylabel('Tail y-position [cm]');
    yt = yticks;
    yl = ylim;
    % --- Right y-axis: scaled to % of body length
    yyaxis right
    ax = gca; % get current axes
    ylim(yl);
    yticks(yt)
    yticklabels(yt /(Ly*100)*100); % set tick labels correctly
    ax.YColor = 'k';
    ylBottomRight = ylabel('\% of body thickness', 'color', 'k');

    % Align left y-labels (left axes)
    if ylTop.Position(1) < ylBottom.Position(1)
        set(ylBottom, 'Position', [ylTop.Position(1), ylBottom.Position(2), ylBottom.Position(3)]);
    else
        set(ylTop, 'Position', [ylBottom.Position(1), ylTop.Position(2), ylTop.Position(3)]);
    end

    % Align right y-labels (right axes)
    if ylTopRight.Position(1) < ylBottomRight.Position(1)
        set(ylTopRight, 'Position', [ylBottomRight.Position(1), ylTopRight.Position(2), ylTopRight.Position(3)]);
    else
        set(ylBottomRight, 'Position', [ylTopRight.Position(1), ylBottomRight.Position(2), ylBottomRight.Position(3)]);
    end

    % add rectangle
    annotation('textbox', [0.875 0.839 0.04 0.07],...
        'BackgroundColor', 'none', ...
        'EdgeColor', 'black',...
        'LineWidth', 1.5)
end