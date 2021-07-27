function plotTriangulation ( c4n, n4e )
%% Draw a triangular grid into a new figure.
%   plotTriangulation(c4n, n4e) draws the grid defined by c4n and n4e into
%                               a figure.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    % Set titles of plot and window.
    title({'Mesh plot'; [num2str(size(c4n,1)), ' nodes']});
    % This can be done with triplot but patch is _much_ faster.
    % Get the coordinates for each node of each triangle.
    X1 = c4n(n4e(:,1),1);
    Y1 = c4n(n4e(:,1),2);
    X2 = c4n(n4e(:,2),1);
    Y2 = c4n(n4e(:,2),2);
    X3 = c4n(n4e(:,3),1);
    Y3 = c4n(n4e(:,3),2);
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    % Set the colour each triangle is filled with.
    C = 'white';
    % Draw everything, make sides blue (looks more like triplot).
    patch(X,Y,C,'EdgeColor','black');
    axis equal tight;
    drawnow;
end
