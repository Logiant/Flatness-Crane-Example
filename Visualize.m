function Visualize(p, theta, t, params, fig_id, SAVE_FIGS)

if nargin==5
    SAVE_FIGS = false;
end


l = params.l;
nt = length(t);
for i = 1:params.di:nt
    figure(fig_id);
    clf; hold on;

    
    %draw the rail
    plot([params.p_min, params.p_max], [-.2, -.2], '-k', 'linewidth', 4)
    
    %plot trolley and bob histories
    plot(p(1:i), zeros([1,i]), '--k');
    plot(p(1:i) + l.*sin(theta(1:i)), -l.*cos(theta(1:i)), '--k');
    
    %draw the trolley
    plot(p(i), 0, 'sk', 'markersize', 20, 'markerfacecolor', 'r');
    
    %draw the rope
    plot( [p(i), p(i) + l*sin(theta(i))], [0, -l*cos(theta(i))], ...
        'linewidth', 3, 'color', [0.4, 0.7, 0.3]  )

    
    %draw the bob
    plot(p(i) + l*sin(theta(i)), -l*cos(theta(i)), 'ko', ...
        'markersize', 20, 'markerfacecolor', 'k');

    %axis([-24, 24, -12, 12]);
    %axis equal;
    %text(0, 1, [ num2str(t(i)) , 'seconds'])
    
    axis equal;
    
    ax = gca;
    text(ax.XLim(1)+1, ax.YLim(2)-1, ['t = ', num2str(t(i))])
    
    pause(0.01)
    if SAVE_FIGS
        saveas(gcf, ['figs/crane_', num2str(i, '%04g'), '.png']);
    end
end
end

