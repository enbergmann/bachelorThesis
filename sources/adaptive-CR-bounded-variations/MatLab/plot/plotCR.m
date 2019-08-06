function plotCR(c4n, n4e, x, OPTtitle)
%% Draw a Crouzeix-Raviart-function.
%   plotCR(c4n, n4e, x, OPTtitle) draws the CR-function defined by the grid
%                                 (c4n, n4e) and the basis coefficients (x)
%                                 The input argument OPTtitle is optional,
%                                 it sets the title of the figure. The 
%                                 default value is empty.

    if nargin == 4
            title(OPTtitle);
        else
            title('');
    end
    
    %% Get coordinates for nodes.
    X1 = c4n(n4e(:,3),1);
    Y1 = c4n(n4e(:,3),2);
    X2 = c4n(n4e(:,1),1);
    Y2 = c4n(n4e(:,1),2);
    X3 = c4n(n4e(:,2),1);
    Y3 = c4n(n4e(:,2),2);
    
    %% Translate values for degrees of freedom into values for nodes.
    s4e = computeS4e(n4e);
    W = x(s4e)';    % Get x for each side of each element.
    Z = (ones(3)-2*eye(3))*W;
    
    %% Assemble parameters for the patch function.
    X = [X1'; X2'; X3'];
    Y = [Y1'; Y2'; Y3'];
    % The colour of a triangle is determined by its midpoint.
    C = sum(Z,1)/3;
    
    %% Put everything together.
    % For large numbers of elements, omit the black boundary around each
    % triangle.
    if( size(n4e,1) > 2000 )
        patch(X,Y,Z,C,'EdgeColor','none');
    else
        patch(X,Y,Z,C);
    end
    view(-37.5,30);
    grid on
end