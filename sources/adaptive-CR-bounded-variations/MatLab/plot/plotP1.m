function plotP1(c4n, n4e, x, OPTtitle)
%% Draw a P1-function.
%   plotP1(c4n, n4e, x, OPTtitle) draws the P1-function defined by the grid
%                                 (c4n, n4e) and the basis coefficients (x).
%                                 The input argument OPTtitle is
%                                 optional, it sets the title of the figure.
%                                 The default value is empty.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    if(size(n4e,1)>2000)
        trisurf(n4e,c4n(:,1),c4n(:,2),x,'EdgeColor','none');
    else
        trisurf(n4e,c4n(:,1),c4n(:,2),x);
    end

    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end
    drawnow;
end
