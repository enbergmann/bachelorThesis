function ax = plotP1image(c4n, u, hmin)
%%PLOTP1IMAGE draws a  plot of an S1 function in current figure. It employs
%linear interpolation between nodal values. This functions works for
%rectangular domains only at the moment. 
%   ax = plotP1image(c4n, u, hmin) gets the coordinates c4n, the 
%   coefficients u of the discrete function, and the minimal mesh-size hmin
%   in the triangulation
%
%This function uses MATLAB PDE Toolbox and is *not* compatible with octave!

% Copyright 2021 Philipp Bringmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

    
    %% PARAMETER
    N = fix(1/hmin);

    %% CREATE MESHGRID
    xmin = min(c4n(:,1));
    xmax = max(c4n(:,1));
    ymin = min(c4n(:,2));
    ymax = max(c4n(:,2));

    [X, Y] = meshgrid(linspace(xmin, xmax, N), linspace(ymin, ymax, N));

    %% INTERPOLATE P1 FUNCTION
    F = scatteredInterpolant(c4n(:,1), c4n(:,2), u, 'linear');
    V = F(X, Y);

    %% SURFACE PLOT
    % Create axes object
    ax = gca;
    % Plot
    if numel(V) > 1
        surf(ax, X, Y, V, 'EdgeColor', 'none');
    else
        return
    end

    %% FORMAT PLOT
    title(ax, {'Interpolated P1 image'; [num2str(size(c4n,1)), ' nodes']});
    colormap(ax, 'gray');
    colorbar(ax);
    axis('equal', 'tight');
    view(ax, 2);
    drawnow;
end

