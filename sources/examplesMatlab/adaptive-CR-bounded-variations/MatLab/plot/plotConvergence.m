function plotConvergence(nrDoF4lvl, error4lvl, OPTname, OPTlimits)
%% Plot error for different levels.
%   plotConvergence(nrDoF4lvl, error4lvl) plots the error given by
%           error4lvl over the degrees of freedom given by nrDoF4lvl, both
%           using logarithmic scale. The optional String argument together
%			with the convergence rate is added to the legend.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
    if nargin > 2
        name = OPTname;
    else
        name = '';
    end
    if nargin > 3
        limits = OPTlimits;
    else
        limits = [1:length(nrDoF4lvl)];
    end

    %% Plot.
    plot = loglog(nrDoF4lvl,error4lvl,'-s');
    title('Convergence Plot');

    %% Set name and title of the figure.
    % If there is only one set of data points, compute the average slope of
    % the curve.
    if isvector(error4lvl)
        if(length(error4lvl)>1)
            p = polyfit(log(nrDoF4lvl(limits)),log(error4lvl(limits)),1);
            % Set convergence rate as title.
            name = char([name ' (' num2str(-p(1)) ')']);
            set(plot,'DisplayName',name);
        end
    else
        set(plot,'DisplayName',name);
    end

    legend('off');
    legend('show');

    drawnow;
end
