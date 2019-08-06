function plotP3(c4n, n4e, x, OPTtitle)
%% Draw a P1-function.
%   plotP3(c4n, n4e, x, OPTtitle) draws the P1-function defined by the grid
%                                 (c4n, n4e) and the basis coefficients (x).
%                                 The input argument OPTtitle is
%                                 optional, it sets the title of the figure.
%                                 The default value is empty.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

  x = x(1:size(c4n,1));
  if size(x, 2) == 1
    if(size(n4e,1)>2000)
      trisurf(n4e,c4n(:,1),c4n(:,2),x,'EdgeColor','none');
    else
      trisurf(n4e,c4n(:,1),c4n(:,2),x);
    end
  elseif size(x, 2) == 2
    set(gca, 'XLim', 1.1*[min(c4n(:,1)), max(c4n(:,1))],...
             'YLim', 1.1*[min(c4n(:,2)), max(c4n(:,2))]);
    quiver2(c4n(:,1), c4n(:,2), x(:,1), x(:,2), 'n=', 0.1, 'w=', [1 1]);
  else
    error('Invalid shape of coefficient vector');
  end

  if nargin == 4
    title(OPTtitle);
  else
    title('');
  end

  drawnow;
end
