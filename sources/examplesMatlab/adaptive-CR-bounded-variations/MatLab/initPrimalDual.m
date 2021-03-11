function [STIMACR, MASSCR, intFSq4e, intFCR14e, intFCR24e, intFCR34e,...
            gradsCR4e, s4e, dof, nrSides, area4e] = initPrimalDual(c4n, n4e, n4sDb, f, degreeF)

    nrElems = size(n4e,1);
    area4e = computeArea4e(c4n,n4e);
    s4e = computeS4e(n4e);
    nrSides = max(max(s4e));
    s4n = computeS4n(n4e);
    
    DbSides = zeros(1,size(n4sDb,1));
    for i = 1:size(n4sDb,1)
        DbSides(i) = s4n(n4sDb(i,1),n4sDb(i,2));
    end
    
    dof = setdiff(1:nrSides,DbSides);
    
    STIMACR4e = zeros(3, 3, nrElems);
    MASSCR4e = zeros(3, 3, nrElems);
    gradsCR4e = zeros(3, 2, nrElems);
    
    for elem = 1 : nrElems
        nodes = n4e(elem,:);   % nodes of this element
        coords = c4n(nodes,:); % coordinates for the nodes
        area = area4e(elem);   % area of this element
        gradsNC = [1 1 1; coords']\[0 0; -2 0; 0 -2]; % gradients for CR basis
        gradsNC = gradsNC([3 1 2],:); % reorder to fit DoF numbering
        gradsCR4e(:,:,elem) = gradsNC;
        STIMACR4e(:,:,elem) = area*(gradsNC*gradsNC'); % local stiffness matrix
        MASSCR4e(:,:,elem) = area*eye(3)/3;
    end
    
    s4eT = s4e';
    I = [s4eT;s4eT;s4eT];
    J = [s4eT(:),s4eT(:),s4eT(:)]';
    
    STIMACR = sparse(I(:),J(:),STIMACR4e(:));
    MASSCR = sparse(I(:),J(:),MASSCR4e(:));
    
    [Lambda1, Lambda2, Lambda3] = barycentricCoords(c4n,n4e);
    
    intFCR14e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)...
        (f(Gpts4p).*(1 - 2*Lambda3(Gpts4p))), degreeF + 1, area4e);
    intFCR24e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)...
        (f(Gpts4p).*(1 - 2*Lambda1(Gpts4p))), degreeF + 1, area4e);
    intFCR34e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)...
        (f(Gpts4p).*(1 - 2*Lambda2(Gpts4p))), degreeF + 1, area4e);
    intFSq4e = integrate(c4n, n4e, @(n4p,Gpts4p,Gpts4ref)f(Gpts4p).^2, 2*degreeF, area4e);
    
end

function [Lambda1,Lambda2,Lambda3]=barycentricCoords(c4n,n4e)
%% barycentricCoords
% This subfunction computes the barycentric coordinates of the given
% triangulation.

    %% Compute a0,a1,..,a4
    % vertices
    a1=c4n(n4e(:,1),:);
    a2=c4n(n4e(:,2),:);
    a3=c4n(n4e(:,3),:);
    %% Define barycentric coordinates
    % by means of affine transformation Phi x + a3
    % where Phi=[a1-a3,a2-a3] '2x2-matrix' in row notation
    % where [a c; b d] --> [a b c d]
    Phi=[a1-a3,a2-a3];
    % Compute Psi=Phi_inverse
    Psi=repmat(1./(Phi(:,1).*Phi(:,4) - Phi(:,2).*Phi(:,3)),1,4)...
          .*[Phi(:,4),-Phi(:,2),-Phi(:,3),Phi(:,1)];
    % Barycentric coordinates
    Lambda1=@(z)(Psi(:,1).*(z(:,1)-a3(:,1)) + Psi(:,3).*...
                                            (z(:,2)-a3(:,2)));
    Lambda2=@(z)(Psi(:,2).*(z(:,1)-a3(:,1)) + Psi(:,4).*...
                                            (z(:,2)-a3(:,2)));
    Lambda3=@(z)(1 - Lambda1(z) - Lambda2(z));
end