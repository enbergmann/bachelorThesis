function plotP04s(c4n, n4e, x, OPTtitle)
%% Draw the P0-function on sides.
%   plotP04s(c4n, n4s, x, OPTtitle) draws the P0-function defined
%                                   by the grid (c4n, n4s) and the coeffi-
%                                   cient vector x into the active figure and
%                                   The input argument OPTtitle is optional, it
%                                   sets the title of the figure. The
%                                   default value is empty.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    angle = [-37.5, 30];  % set the point of view


    n4s = computeN4s(n4e);
    [errX,errY,errZ,errC] = getErrorXYZC(c4n, n4s, x);
    if(size(n4e,1)>2000)
    	patch(errX,errY,errZ,errC,'EdgeColor','none');
    else
    	patch(errX,errY,errZ,errC);
    end

    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end

    view(angle);
    grid on
    drawnow;
end

%% Function to get the coordinates for the patch-function
function [valX, valY, valZ,valC] = getErrorXYZC(c4n, n4s, x)
    % coordinates for all nodes of sides
    X1 = c4n(n4s(:,1), 1);
    X2 = c4n(n4s(:,2), 1);
    Y1 = c4n(n4s(:,1), 2);
    Y2 = c4n(n4s(:,2), 2);

    % sides for the function in patch-style
    valX = [X1';X1';X2';X2'];
    valY = [Y1';Y1';Y2';Y2'];

    % each side should have the same height - the function value
    valZ = [zeros(size(x'));x';x';zeros(size(x'))];
    valC = valZ / max(x);
end
