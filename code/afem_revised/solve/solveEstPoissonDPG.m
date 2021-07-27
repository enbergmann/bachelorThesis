function [ u, ndof, u4e, globalErr, localErr ] = solveEstPoissonDPG(f, u4Db, c4n, n4e, n4sDb, degreeF,...
    baseU, gradU, nCPu, nIPu, c4BaseU, baseQ, baseV, gradV, nCPv, nIPv, kU, kQ, kV)
%SOLVEPOISSONDPG Solves Poisson model problem with the primal dPG method.
%    u = SOLVEESTPOISSONDPG(f, u4Db, c4n, n4e, n4sDb, degreeF, baseU, gradU,
%    nEPu, nIPu, c4BaseU, baseQ, baseV, gradV, nEPv, nIPv, kU, kQ, kV) returns 
%    the DPG coefficient vector of the DPG solution to the Poisson model
%    problem on the triangulation c4n and n4e for given right-hand side f
%    and boundary data u4Db on the Dirichlet boundary n4sDb and given polynomial
%    basis for the test space and trial space. degreeF is the
%    polynomial degree used for the Gaussian quadrature rule.
%
% input:    f           - right-hand side           
%           u4Db        - Dirichlet boundary data
%           c4n         - coordinates for nodes
%           n4e         - nodes for element
%           n4sDb       - nodes for Dirichlet boundary sides
%           degreeF     - polynomial degree for Gaussian quadrature rule
%                         for f
%           baseU       - basis for P_{k_u}(\mathcal{T}) as a part of the
%                         trial space, for simplicity only on the reference
%                         triangle; order: base functions in vertices,
%                         on first side, second side, third side, inner
%                         points
%           gradU       - gradients of baseU (in the same order)
%           nCPu        - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIPu        - number of base functions in inner points
%           c4BaseU     - array of values which determine the position of
%                         the base functions on the side, e.g. for P2 [0.5] 
%           baseV       - basis for P_{k_v}(\mathcal{T}) as the test space,
%                         for simplicity only on the reference triangle;
%                         same order as for baseU
%           baseQ       - basis for P_{k_v}(\mathcal{E}) as a part of the
%                         trial space, for simplicity only on the reference
%                         side
%           gradV       - gradients of baseV (in the same order)
%           nCPv        - 1 if there are base functions in vertices, 0
%                         otherwise
%           nIPv        - number of base functions in inner points
%           kU          - total polynomial degree of P_{k_u}(\mathcal{T})
%           kQ          - polynomial degree for P_{k_q}(\mathcal{E}) as
%                         a part of the trial space
%           kV          - total polynomial degree of P_{k_v}(\mathcal{T})
%
% output:   u           - solution to the Primal DPG method for the 
%                         Poisson model problem, order: coefficients for 
%                         baseU, coefficients for baseQ 
%           ndof        - number of degrees of freedom
%           u4e         - solution to the Primal DPG method for the Poisson
%                         model problem for every triangle  
%           globalErr   - estimated error
%           localErr    - estimated error on each triangle
%


%% INITIALIZATION
n4s = computeN4s(n4e);
s4e = computeS4e(n4e);
s4n = computeS4n(n4e);
nElem = size(n4e, 1);
nNodes = size(c4n, 1);
nSides = size(n4s, 1);

nBaseU = size(baseU, 2);
nBaseQ = kQ+1;
nBaseV = size(baseV,2);

nSPu = (nBaseU - 3*nCPu - nIPu)/3;
nSPv = (nBaseV - 3*nCPv - nIPv)/3;

dimU = nNodes*nCPu+nSides*nSPu+nElem*nIPu;
dimV = nElem*nBaseV;
dimQ = nSides*nBaseQ;

area4e = computeArea4e(c4n, n4e);
length4s = computeLength4s(c4n,n4s);

signs4e = computeSigns4e(c4n,n4e, n4s, s4e);
indicatorPos4e = (1+signs4e)./2; % 1 if positive sign in signs4e else 0
indicatorNeg4e = (1-signs4e)./2; % 1 if negative sign in signs4e else 0

%% DEGREES OF FREEDOM
DbNodes = unique(n4sDb);
DbSides = full(s4n(sub2ind(size(s4n), n4sDb(:,1), n4sDb(:,2))));
dof = 1:(dimU+dimQ);

if nCPu == 1    
    dof = setdiff(dof,DbNodes);
end

if nSPu ~=0    
   dof = setdiff(dof, nNodes*nCPu + repmat(nSPu*(DbSides-1),1,nSPu) + ...
            repmat(1:nSPu,size(DbSides,1),1));    
end

ndof = length(dof);

%% GET LOCAL SHAPE FUNCTIONS AND GRADIENTS
[uPhi4e, uGradPhi4e] = transformAffine(c4n, n4e, baseU, gradU);
[vPhi4e, vGradPhi4e] = transformAffine(c4n, n4e, baseV, gradV);

%% BUILD LOCAL STIFFNESS MATRIX

% local stiffness matrix for integral(grad(u)*grad(v)) on each triangle
uB4e = nan(nBaseV, nBaseU, nElem);
for j = 1:nBaseV
  for k = 1:nBaseU
    % evaluate bilinear form
    INT = @(n4p, Gpts4p, Gpts4ref)...
           sum(uGradPhi4e{k}(Gpts4ref).*vGradPhi4e{j}(Gpts4ref), 2);
    uB4e(j,k,:) = integrate(INT, c4n, n4e, kU+kV-2, area4e);
  end
end

% local stiffness matrix for integral(t*v) on the trace of each triangle
qB4e = nan( nElem, 3*nBaseQ, nBaseV);
for j=1:nBaseV
    for k = 1:nBaseQ    
        % for each side of the triangle
        INT1 = @(n4p, Gpts4p, Gpts4ref)...
            baseQ{k}(Gpts4ref).*baseV{j}([zeros(size(Gpts4ref)),1-Gpts4ref]);
        INT2 = @(n4p, Gpts4p, Gpts4ref)...
            baseQ{k}(Gpts4ref).*baseV{j}([Gpts4ref, zeros(size(Gpts4ref))]);
        INT3 = @(n4p, Gpts4p, Gpts4ref)...
            baseQ{k}(Gpts4ref).*baseV{j}(kron(Gpts4ref,[-1, 1]) + ... 
            repmat([1 0], size(Gpts4ref)));
        
        temp1 = integrate(INT1, c4n, n4e(:,[2 3]), kQ+kV, length4s(s4e(:,1)));
        temp2 = integrate(INT2, c4n, n4e(:,[3 1]), kQ+kV, length4s(s4e(:,2)));
        temp3 = integrate(INT3, c4n, n4e(:,[1 2]), kQ+kV, length4s(s4e(:,3)));
        temp = signs4e .* [temp1 temp2 temp3];
        
        qB4e(:,k+[0 nBaseQ 2*nBaseQ],j) = -temp;
    end
end
qB4e = permute(qB4e , [3 2 1]);

%% BUILD LOCAL NORM MATRIX
M4e = nan(nElem, nBaseV, nBaseV);
for j = 1:nBaseV
    for k= j:nBaseV
        INT = @(n4p, Gpts4p, Gpts4ref)...
           baseV{j}(Gpts4ref).*baseV{k}(Gpts4ref);
        M4e(:,j,k) = integrate(INT,c4n,n4e, 2*kV, area4e);
        INT = @(n4p, Gpts4p, Gpts4ref)...
           sum(vGradPhi4e{j}(Gpts4ref).*vGradPhi4e{k}(Gpts4ref), 2);        
        M4e(:,j,k) = M4e(:,j,k)+integrate(INT,c4n,n4e, 2*kV-2, area4e);
        if k > j
            M4e(:,k,j) = M4e(:,j,k);  % use symmetry
        end
    end
end
M4e = permute(M4e, [3 2 1]);

%% BUILD LOCAL RIGHT-HAND SIDE
res = cell(1,nBaseV);
for i=1:nBaseV
    res{i} = fun(baseV{i}, c4n, n4e, degreeF, area4e, kV, f);    
end
F4e = cell2mat(res);

%% ASSEMBLING RIGHT HAND SIDE
F = F4e';
F = F(:);

%% CONVERT LOCAL TO GLOBAL DOFS FOR NORM MATRIX
invM4e = nan(nBaseV, nBaseV, nElem);
Idof4e = zeros(nBaseV, nBaseV, nElem);
for k=1:nElem
    invM4e(:,:,k) = M4e(:,:,k)\eye(nBaseV);
    Idof4e(:,:,k) = repmat((1:nBaseV)+(k-1)*nBaseV, nBaseV,1);
end
Jdof4e = permute(Idof4e, [2 1 3]);
invM = sparse(Idof4e(:), Jdof4e(:), invM4e(:));

%% CONVERT LOCAL TO GLOBAL DOFS FOR STIFFNESS MATRIX
uIdof4e = nan(nBaseV, nBaseU, nElem);
uJdof4e = nan(nBaseV, nBaseU, nElem);

qIdof4e = nan(nBaseV, 3*nBaseQ, nElem);
qJdof4e = nan(nBaseV, 3*nBaseQ, nElem);

%   o --------- o   
%   | \         |   For more than 1 base function on each side, the order 
%   |   a   2   |   in the local stiffness matrix is [b a] on triangle 1  
%   |     \     |   and[a b] on 2. Therefore mapU and mapQ compute the
%   |   1   b   |   right order for the respective basis.
%   |         \ |
%   o --------- o    

mapU = nCPu*nNodes + kron(nSPu*(s4e(:) - indicatorPos4e(:)) + ...
    indicatorNeg4e(:), ones(1,nSPu)) +  kron(signs4e(:),1:nSPu);
mapQ = dimU + kron(nBaseQ*(s4e(:)-indicatorPos4e(:)) + ...
    indicatorNeg4e(:), ones(1,nBaseQ)) + kron(signs4e(:),1:nBaseQ);

for k = 1:nElem
    uIdof4e(:,:,k) = repmat((1:nBaseV)'+ (k-1)*nBaseV,1,nBaseU);
    
    temp = nan(1,nBaseU);
    if nCPu == 1
       temp([1 2 3]) = n4e(k,:);        
    end
    
    if nSPu ~=0
        temp(nCPu*3+(1:(3*nSPu))) = [mapU(k,:), mapU(nElem+k,:), ...
            mapU(2*nElem+k,:)];
    end
    
    if nIPu ~= 0
       temp((nCPu*3+3*nSPu+1):end) = nCPu*nNodes+nSPu*nSides + ...
           (k-1)*nIPu+(1:nIPu);  
    end
    
    uJdof4e(:,:,k) = repmat(temp, nBaseV,1);
    
    qIdof4e(:,:,k) = repmat((1:nBaseV)' + (k-1)*nBaseV,1,3*nBaseQ);    
    qJdof4e(:,:,k) = repmat([mapQ(k,:), mapQ(k+nElem,:), ...
        mapQ(k+2*nElem,:)], nBaseV, 1);
end

%% ASSEMBLING OF STIFFNESS MATRIX
Idof = [uIdof4e(:); qIdof4e(:)];
Jdof = [uJdof4e(:); qJdof4e(:)];
B = sparse(Idof, Jdof, [uB4e(:); qB4e(:)], dimV, dimU + dimQ);

A = B'*invM*B; 
C = B'*invM*F;

%% INCLUDE DIRICHLET BOUNDARY CONDITIONS
u = zeros(dimU+dimQ,1);

c4nDbT = c4n(DbNodes,:);

if nCPu==1
    u(DbNodes) = u4Db(c4nDbT);
end

if nSPu~=0
    nrDbSides = size(DbSides,1);
    
    % coordinates for the base functions on the border sides
    temp1 = kron(c4n(n4s(DbSides,1),:),1-c4BaseU); 
    temp2 = kron(c4n(n4s(DbSides,2),:),c4BaseU);
    c4nDbE = reshape(temp1,nSPu*nrDbSides,2)+reshape(temp2,nSPu*nrDbSides,2); 
    
    u(nNodes*nCPu + reshape(repmat(nSPu*(DbSides-1),1,nSPu) + ...
              repmat(1:nSPu,nrDbSides,1),nSPu*nrDbSides,1)) = u4Db(c4nDbE);
end

C = C - A*u;

%% SOLVE LINEAR SYSTEM
u(dof) = A(dof,dof) \ C(dof);

%% COMPUTE LOCAL COEFFICIENTS
u4e = u(uJdof4e(1,:,:));
u4e = permute(u4e, [3,2,1]);

%% COMPUTE LOCAL ERROR
temp = F-B*u;
localErr = nan(nElem,1);
for k = 1:nElem
    localErr(k) = temp((k-1)*nBaseV+(1:nBaseV))' * invM4e(:,:,k) * ...
                  temp((k-1)*nBaseV+(1:nBaseV));
end

%% COMPUTE GLOBAL ERROR
globalErr = temp' * invM*temp; % sum(localErr)

localErr = sqrt(localErr);
globalErr = sqrt(globalErr);
end

function X =fun(Phi, c4n, n4e, degreeF, area4e, kV,f)
    INT = @(n4p,Gpts4p,Gpts4ref) f(Gpts4p).*Phi(Gpts4ref);    
    X = integrate(INT, c4n, n4e, degreeF+kV, area4e);
end