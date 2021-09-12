function plotInitialTriangulation ( c4n, n4e )
%% Draw a triangular grid into a new figure.
%   plotTriangulation(c4n, n4e) draws the grid defined by c4n and n4e into
%                               a figure.

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
    patch(X,Y,C,'EdgeColor','blue');
    
    nrNodes = size(c4n,1);
    nrs = num2str( (1:nrNodes)');
    text(c4n(:,1),c4n(:,2),nrs,...
            'HorizontalAlignment','Center','BackgroundColor','White',...
            'EdgeColor','Blue','FontSize',16);
    
    mid4e = computeMid4e(c4n,n4e);
    nrElems = size(n4e,1);
    latexNrs = [repmat('$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {',...
                nrElems,1),num2str( (1:nrElems)'),repmat('}}}$',nrElems,1)];
    text(mid4e(:,1),mid4e(:,2),latexNrs,'Interpreter','latex',...
      'Color','Blue','FontSize',18);

    n4s = computeN4s(n4e);
    mid4s = computeMid4s(c4n,n4s);
    ns = size(mid4s,1);
    latexNrs = [repmat('$\raisebox{.5pt}{\raisebox{-.9pt} {',...
                ns,1),num2str((1:ns)'),repmat('}}$',ns,1)];
    text(mid4s(:,1),mid4s(:,2),latexNrs,'Interpreter','latex',...
      'Color','Blue','FontSize',18);
    drawnow;
end