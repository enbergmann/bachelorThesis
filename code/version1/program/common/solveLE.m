function u = solveLE(c4n,n4e,n4sDb,n4sNb,f,Lambda,alpha,)
  % return STIMANC: Steifigkeitsmatrix
  % return MAMANC: Massematrix

    %% Initialisation
    nrElems = size(n4e,1);
    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));
    area4e = computeArea4e(c4n,n4e);
    mid4e = computeMid4e(c4n,n4e);
    s4n = computeS4n(n4e);
    % Dirichlet boundary sides
    DbSides = zeros(1,size(n4sDb,1));
    for i = 1:size(n4sDb,1)
        DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
    end
    % Neumann boundary sides
    NbSides = zeros(1,size(n4sNb,1));
    for i = 1:size(n4sNb,1)
        NbSides(i) = s4n(n4sNb(i,1),n4sNb(i,2));
    end
    % degrees of freedom: one per non-Dirichlet side
    dof = setdiff(1:nrSides,DbSides);
    nrDof = length(dof);
    
    STIMANClocal = zeros(3,3,nrElems);
    MAMANClocal = repmat(eye(3)/3,1,1,nrElems);
    b = zeros(nrSides,1);

    %% Create the stiffness matrix, mass matrix and right-hand side b
    for elem = 1 : nrElems
        nodes = n4e(elem,:);   % nodes of this element
        sides = s4e(elem,:);   % sides of this element
        coords = c4n(nodes,:); % coordinates for the nodes
        area = area4e(elem);   % area of this element
        grads = [1 1 1; coords']\[0 0; -2 0; 0 -2]; % gradients for CR basis
        grads = grads([3 1 2],:); % reorder to fit DoF numbering
        STIMANClocal(:,:,elem) = area * (grads * grads'); % local stiffness matrix
        MAMANClocal(:,:,elem) = area * MAMANClocal(:,:,elem); % local mass matrix
        mid = mid4e(elem,:);     % midpoint of this element
     %   b(sides) = b(sides) + area*f(mid)*ones(3,1)/3; % right-hand side
        INT = @(n4p,Gpts4p,Gpts4ref) f(Gpts4p).*;
        temp =
        b(sides) = b(sides) + integrate(INT,c4n,n4e,6,area4e)(sides);
    end

    % assembly of the global stiffness matrix and global mass matrix
    s4eT = s4e';
    I = [s4eT;s4eT;s4eT];
    J = [s4eT(:),s4eT(:),s4eT(:)]';
    STIMANC = sparse(I(:),J(:),STIMANClocal(:));
    MAMANC = sparse(I(:),J(:),MAMANClocal(:));
    
    %% Solve System
    A = STIMANC/tau+alpha*MAMANC; 
end