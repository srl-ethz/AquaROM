% ------------------------------------------------------------------------ 
% fig_head_tail_motion.m
%
% Description: create a figure showing the motion of the head and the tail
% node.
%
% Last modified: 08/06/2025, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function fig = fig_head_tail_motion(uHead_sol, uTail_sol, timePlot, ...
    Lx, Ly, rectPos)

    fig = figure('units','centimeters','position',[3 3 9 7.0]);

    % x-position (Head) ___________________________________________________
    uHead_FOM = uHead_sol.uHead_FOM(1,:);
    uHead_ROM = uHead_sol.uHead_ROM(1,:);
    uHead_PROM_3 = uHead_sol.uHead_PROM_3(1,:);
    uHead_PROM_5 = uHead_sol.uHead_PROM_5(1,:);
    uHead_PROM_8 = uHead_sol.uHead_PROM_8(1,:);
    
    subplot(2,1,1);
    hold on;
    p1 = plot(timePlot, uHead_FOM, '--', 'DisplayName', 'FOM', 'LineWidth', 1.0);
    p2 = plot(timePlot, uHead_ROM, 'DisplayName', 'ROM', 'LineWidth', 1.0);
    p3 = plot(timePlot, uHead_PROM_3, '-.', 'DisplayName', 'PROM 3p', 'LineWidth', 1.0);
    p4 = plot(timePlot, uHead_PROM_5, '-.', 'DisplayName', 'PROM 5p', 'LineWidth', 1.0);
    p5 = plot(timePlot, uHead_PROM_8, '-.', 'DisplayName', 'PROM 8p', 'LineWidth', 1.0);
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
    yticklabels(yt /(Lx*100)*100); % set tick labels correctly
    ax.YColor = 'k';
    ylTopRight = ylabel('\% of body length', 'color', 'k');

    % --- Back to left for legend
    yyaxis left
    lgd = legend([p1,p2,p3,p4,p5],{'FOM', 'ROM',  'PROM 3p', 'PROM 5p', 'PROM 8p'},'Location','northwest', 'NumColumns',1, 'interpreter', 'latex');
    lgd.Position(2) = lgd.Position(2) + 0.08;

    % y-position (Tail) ___________________________________________________
    uTail_FOM = uTail_sol.uTail_FOM(2,:);
    uTail_ROM = uTail_sol.uTail_ROM(2,:);
    uTail_PROM_3 = uTail_sol.uTail_PROM_3(2,:);
    uTail_PROM_5 = uTail_sol.uTail_PROM_5(2,:);
    uTail_PROM_8 = uTail_sol.uTail_PROM_8(2,:);

    subplot(2,1,2);
    hold on;
    plot(timePlot, uTail_FOM, '--', 'DisplayName', 'FOM', 'LineWidth', 1.0);
    plot(timePlot, uTail_ROM, 'DisplayName', 'ROM', 'LineWidth', 1.0);
    plot(timePlot, uTail_PROM_3, '-.', 'DisplayName', 'PROM 3p', 'LineWidth', 1.0);
    plot(timePlot, uTail_PROM_5, '-.', 'DisplayName', 'PROM 5p', 'LineWidth', 1.0);
    plot(timePlot, uTail_PROM_8, '-.', 'DisplayName', 'PROM 8p', 'LineWidth', 1.0);
    grid on;
    xlabel('Time [s]');
    
    % --- Left y-axis: raw cm
    yyaxis left
    ylBottom = ylabel('Tail y-position [cm]');
    yt = yticks;
    yl = ylim;

    % --- Right y-axis: scaled to % of body thickness
    yyaxis right
    ax = gca; % get current axes
    ylim(yl);
    yticks(yt)
    yticklabels(yt /(Ly*100)*100); % set tick labels correctly
    ax.YColor = 'k';
    ylBottomRight = ylabel('\% of body thickness', 'color', 'k');

    % Axis Alignment ______________________________________________________
    % left y-labels (left axes)
    if ylTop.Position(1) < ylBottom.Position(1)
        set(ylBottom, 'Position', [ylTop.Position(1), ylBottom.Position(2), ylBottom.Position(3)]);
    else
        set(ylTop, 'Position', [ylBottom.Position(1), ylTop.Position(2), ylTop.Position(3)]);
    end

    % right y-labels (right axes)
    if ylTopRight.Position(1) < ylBottomRight.Position(1)
        set(ylTopRight, 'Position', [ylBottomRight.Position(1), ylTopRight.Position(2), ylTopRight.Position(3)]);
    else
        set(ylBottomRight, 'Position', [ylTopRight.Position(1), ylBottomRight.Position(2), ylBottomRight.Position(3)]);
    end

    % add rectangle
    if ~isempty(rectPos)
        annotation('textbox', rectPos,...
        'BackgroundColor', 'none', ...
        'EdgeColor', 'black',...
        'LineWidth', 1.5)
    end
    
end
