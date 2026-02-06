% ------------------------------------------------------------------------ 
% fig_head_tail_motion.m
%
% Description: create a figure showing the motion of the head and the tail
% node.
%
% Last modified: 06/02/2026, Mathieu Dubied, ETH Zurich
% ------------------------------------------------------------------------
function fig = fig_comparison_PROM_approximation(timePlot, ...
                                uHead_exact, uHead_approx, ...
                                uTail_exact, uTail_approx, ...
                                Lx, Ly, position)

    fig = figure('units','centimeters','position', position);

    % x-position (Head) ___________________________________________________
    uHead_exact = uHead_exact(1,:);
    uHead_approx = uHead_approx(1,:);
    
    subplot(2,1,1);
    hold on;
    p1 = plot(timePlot, uHead_exact, 'DisplayName', 'Ground truth', 'LineWidth', 1.0);
    p2 = plot(timePlot, uHead_approx, '--','DisplayName', 'Approximation', 'LineWidth', 1.0);
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
    lgd = legend('Location','northwest', 'NumColumns',1, 'interpreter', 'latex');
    lgd.Position(2) = lgd.Position(2) + 0.08;

    % y-position (Tail) ___________________________________________________
    uTail_exact = uTail_exact(2,:);
    uTail_approx = uTail_approx(2,:);

    subplot(2,1,2);
    hold on;
    p1 = plot(timePlot, uTail_exact, 'DisplayName', 'Ground truth', 'LineWidth', 1.0);
    p2 = plot(timePlot, uTail_approx, '--','DisplayName', 'Approximation', 'LineWidth', 1.0);
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
    
end
