function [BAGradJN4e, BAGradJN4e2, NormGradJN4e]...
            = computeBAGradJN(u, gradCR4e, m, p, r, t, c4n, n4e, n4sDb, s4e, nrSides, area4e)
%% computeBAGradJN - compute the norm of Jn - Pi0 Jn
%
% Input:     u              solution vector
%            gradCR4e       piecewise gradient of discrete solution
%            m              number of components of f
%            p,r,t          parameters
%            c4n            coordinates for the nodes of the mesh
%            n4e            nodes for the elements of the mesh
%            s4e            sides for the elements
%            n4sDb          nodes for the boundary sides
%            nrSides        number of sides
%            area4e         area of each element
%
% Output:    BAGradJN4e     Lp norm of (1 - Pi0) Grad Jn for each element
%            BAGradJN4e2    L r/(r-t) norm of (1 - Pi0) Grad Jn for each element
%            NormGradJN4e   Lp norm of Grad Jn for each element

    nrElems = size(n4e, 1);
    nrNodes = size(c4n, 1);
    
    gradConf4e = zeros(3,2,nrElems);
    coords4e = zeros(3,2,nrElems);
    
    for j = 1 : nrElems
        nodes = n4e(j,:);   % nodes of the triangle
        coords = c4n(nodes,:); % coordinates of the three nodes
        coords4e(:,:,j) = coords;
        gradConf4e(:,:,j) = [1,1,1;coords'] \ [0,0;eye(2)];
    end
    
    j14n = computeJ14n(u, c4n, n4e, s4e, n4sDb, nrSides, m);
    gradJ14e = gradJ1forE(j14n,gradConf4e,n4e,m,nrNodes);
    
    u = reshape(u,nrSides,m)';
    j14n = reshape(j14n,nrNodes,m)';
    
    uloc4e = zeros(m,3,nrElems);
    j1loc4e = zeros(m,3,nrElems);
    
    for j = 1:nrElems
        uloc4e(:,:,j) = u(:,s4e(j,:));
        j1loc4e(:,:,j) = j14n(:,n4e(j,:));
    end
    
    BAGradJN4e = zeros(nrElems, 1);
    BAGradJN4e2 = zeros(nrElems, 1);
    NormGradJN4e = zeros(nrElems, 1);
    
    for j = 1:nrElems % element-wise computation of the norm of (1 - Pi0) DJn
        coords =  coords4e(:,:,j);
        area = area4e(j);
        gradJ1 = gradJ14e(j,:);
        uloc = uloc4e(:,:,j);
        j14nloc = j1loc4e(:,:,j);
        gradJnUCR = @(x)gradJn(x, uloc, gradJ1, j14nloc, coords, m);
        gradUCR = gradCR4e(j,:);
        
        BAGradJN4e(j) = integrate(coords,[1, 2, 3],...
                 @(n4p,Gpts4p,Gpts4ref)((sum((gradJnUCR(Gpts4p)-gradUCR).^2,2))^(p/2)),10,area);
        NormGradJN4e(j) = integrate(coords,[1, 2, 3],...
                 @(n4p,Gpts4p,Gpts4ref)((sum(gradJnUCR(Gpts4p).^2,2))^(p/2)),10,area);
        if abs(p-r/(r-t)) > 1e-8
            BAGradJN4e2(j) = integrate(coords,[1, 2, 3],...
                 @(n4p,Gpts4p,Gpts4ref)((sum((gradJnUCR(Gpts4p)-gradUCR).^2,2))^(r/(2*(r-t)))),10,area);
        end
    end
    if abs(p-r/(r-t)) < 1e-8
        BAGradJN4e2 = BAGradJN4e;
    end
end

function gradJ14e = gradJ1forE(j14n,gradConf4e,n4e,m,nrNodes)
%% compute the piecewise gradient of J1
    nrElems = size(n4e,1);
    gradJ14e = zeros(nrElems,2*m);
    for j = 1:nrElems
        gradConfloc = gradConf4e(:,:,j);
        gradJ1E = zeros(1,2*m);
        for k = 1:m
            coords = j14n((k-1)*nrNodes + n4e(j,:));
            gradJ1E([2*k-1,2*k]) = coords(1)*gradConfloc(1,:) + coords(2)*gradConfloc(2,:) + coords(3)*gradConfloc(3,:);
        end
        gradJ14e(j,:) = gradJ1E;
    end
end

function val = gradJn(x, uloc, gradJ1, j14nloc, coords, m)
%% gradJn - gradient of Jn for a given element

    nrComps = size(x,1);
    val = zeros(nrComps,2*m);
    grads = [1,1,1;coords'] \ [0,0;eye(2)]; % gradient of nodal basis functions
    
    % compute barycentric coordinates
    a1 = coords(1,:);
    a2 = coords(2,:);
    a3 = coords(3,:);
    Phi = [a1-a3,a2-a3];
    Psi = 1/(Phi(1)*Phi(4) - Phi(2)*Phi(3))*[Phi(4), -Phi(2), -Phi(3), Phi(1)];
    
    Lambda1 = Psi(1)*(x(:,1)-a3(1)) + Psi(3)*(x(:,2)-a3(2));
    Lambda2 = Psi(2)*(x(:,1)-a3(1)) + Psi(4)*(x(:,2)-a3(2));
    Lambda3 = 1 - Lambda1 - Lambda2;
    
    for j = 1:nrComps
        for k = 1:m
            midValues = uloc(k,:) - (j14nloc(k,:) + j14nloc(k,[2,3,1]))/2;
            val(j,[2*k-1,2*k]) = 6*midValues(1)*(Lambda1(j)*grads(2,:) + Lambda2(j)*grads(1,:))...
                                 + 6*midValues(2)*(Lambda2(j)*grads(3,:) + Lambda3(j)*grads(2,:))...
                                 + 6*midValues(3)*(Lambda1(j)*grads(3,:) + Lambda3(j)*grads(1,:));
        end
        val(j,:) = val(j,:) + gradJ1;
    end
end

function j14n = computeJ14n(u, c4n, n4e, s4e, n4sDb, nrSides, m)
%% computeU4n - compute nodal values of J1 of a CR function

    nrNodes = size(c4n,1);
    DbNodes = unique(n4sDb);
    j14n = zeros(1,m*nrNodes);
    nrE4n = zeros(1,nrNodes);
    
    for j = 1:size(n4e,1)
        for k = 1:3
            nrE4n(n4e(j,k)) = nrE4n(n4e(j,k)) + 1;
        end
    end
    
    for j = 1:m
        sides = (j-1)*nrSides + s4e;
        nodes = (j-1)*nrNodes + n4e;
        nodesDb = (j-1)*nrNodes + DbNodes;
        temp = [u(sides(:,1)) + u(sides(:,3)) - u(sides(:,2)); u(sides(:,1))...
            + u(sides(:,2)) - u(sides(:,3)); u(sides(:,2)) + u(sides(:,3)) - u(sides(:,1))]';
        for k = 1:size(n4e,1)
            for ell = 1:3
                j14n(nodes(k,ell)) = j14n(nodes(k,ell)) + temp(k,ell);
            end
        end
        j14n(nodesDb) = zeros(1,length(DbNodes));
        j14n((j-1)*nrNodes + [1:nrNodes]) = j14n((j-1)*nrNodes + [1:nrNodes])./nrE4n;
    end
end

% Copyright 2018 Tran Ngoc Tien
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. No WARRANTY!