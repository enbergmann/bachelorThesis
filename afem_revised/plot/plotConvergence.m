function plotConvergence(nrDoF4lvl, error4lvl, OPTname)
%% Plot error for different levels.
%   plotConvergence(nrDoF4lvl, error4lvl) plots the error given by
%           error4lvl over the degrees of freedom given by nrDoF4lvl, both
%           using logarithmic scale. The optional String argument together
%           with the convergence rate is added to the legend.
    if nargin > 2
        name = OPTname;
    else
        name = '';
    end

    cOrder = [0 0 1; 1 0 0; 0 .7 0; .7 0 .7; .7 .7 0; 0 .7 .7; .3 1 0; 0 .3 1; 0 0 0];

    nlines = get(gcf, 'UserData');

    if ~isempty(nlines)
      nlines = nlines + 1;
    else
      nlines = 1;
    end
    set(gcf, 'UserData', nlines);

    %% Plot.
    plot = loglog(nrDoF4lvl, error4lvl, '-s', 'Color', cOrder(nlines, :));
    title('Convergence history plot');

    %% Set name and title of the figure.
    set(plot,'DisplayName',name);

    legend('off');
    legend('show');
    handle = findobj(gcf,'type','axes','Tag','legend');
    set(handle,'Location','best');

    drawnow;
end

