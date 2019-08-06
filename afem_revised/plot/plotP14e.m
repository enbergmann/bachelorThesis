function plotP14e(c4n, n4e, x, OPTtitle)
%% Draw the P0-function on elements.
%   plotP04e(c4n, n4e, x, OPTtitle) draws the P0-function defined by
%                                   the grid (c4n, n4e) and the coefficient
%                                   vector x into the given figure. The
%                                   input argument OPTtitle is optional, it
%                                   sets the title of the figure. The
%                                   default value is empty.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    %% OPTIONS
    angle = [-37.5, 30];

    %% Get coordinates for nodes.
    X1 = c4n(n4e(:,1),1);
    Y1 = c4n(n4e(:,1),2);
    X2 = c4n(n4e(:,2),1);
    Y2 = c4n(n4e(:,2),2);
    X3 = c4n(n4e(:,3),1);
    Y3 = c4n(n4e(:,3),2);

    %% Translate values for degrees of freedom into values for nodes.
    Z = reshape(x, 3, size(n4e, 1));

    %% Assemble parameters for the patch function.
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    % The colour of a triangle is determined by its midpoint.
    C = sum(Z)/3;

    %% PLOT
    if(size(n4e,1)>2000)
    	patch(X,Y,Z,C,'EdgeColor','none');
    else
    	patch(X,Y,Z,C);
    end

    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end

    view(angle);
    grid on;
    drawnow;
end
