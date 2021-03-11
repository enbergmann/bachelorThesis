function [err,error4e]=error4eStokesLSC(c4n,n4e,f,gradExact1,gradExact2,sigmaExact1,sigmaExact2,uApprox,sigmaApprox)
%% error4eStokesLSC - estimate error of LSFE solution
% using edge-oriented basis
% 
% 
%
% COMMAND: 
%   error4e=error4eStokesRT0Energy(c4n,n4e,gradExact,sigmaApprox)
%
% INPUT:
%   c4n         coordinates for the nodes of the mesh
%   n4e         nodes for the elements of the mesh
%
% OUTPUT:
%   error4e     energy error for each element

    %% Parameter
    % degree for numerical integration
    deg=10;
    
    
    %% Initialization
    % geometry
    n4s=computeN4s(n4e);
    s4e=computeS4e(n4e);
    nrElems=size(n4e,1);
    normal4s=computeNormal4s(c4n,n4s);
    normal4e=computeNormal4e(c4n,n4e);
    area4e=computeArea4e(c4n,n4e);
    
    % functions
    [Lambda1,Lambda2,Lambda3]=barycentricCoords(c4n,n4e);
    
    % integrals
    intF4e=...
        integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p)),deg,area4e);
                            % integral of each component of f 
                            % over this element
    intFSq4e=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)sum(f(Gpts4p).^2,2),deg,area4e);
                            % integral of the scalar product f^2
                            % over this element
                            
    intSigma1Baryc1forE=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact1(Gpts4p).*...
                        [Lambda1(Gpts4p),Lambda1(Gpts4p)]),deg,area4e);
    intSigma1Baryc2forE=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact1(Gpts4p).*...
                        [Lambda2(Gpts4p),Lambda2(Gpts4p)]),deg,area4e);
    intSigma1Baryc3forE=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact1(Gpts4p).*...
                        [Lambda3(Gpts4p),Lambda3(Gpts4p)]),deg,area4e);
    
    intSigma2Baryc1forE=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact2(Gpts4p).*...
                        [Lambda1(Gpts4p),Lambda1(Gpts4p)]),deg,area4e);
    intSigma2Baryc2forE=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact2(Gpts4p).*...
                        [Lambda2(Gpts4p),Lambda2(Gpts4p)]),deg,area4e);
    intSigma2Baryc3forE=...
      integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact2(Gpts4p).*...
                        [Lambda3(Gpts4p),Lambda3(Gpts4p)]),deg,area4e);

    intSigmaSq4e=integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                        sum(sigmaExact1(Gpts4p).*sigmaExact1(Gpts4p) + ...
                        sigmaExact2(Gpts4p).*sigmaExact2(Gpts4p),2),...
                        deg,area4e);

    
    %% Constants
    MA=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],...
        [-4,-2,0,2,4],6,6);
  
   
    %% Variables
    error4e=zeros(nrElems,1);
    contrL2Sigma4e=zeros(nrElems,1);  % L2-norm of sigma - sigma_LS
    contrDivSigma4e=zeros(nrElems,1); % L2-norm of div(sigma - sigma_LS)
    
    
    %% Contributions of the L2-norm of Du - Du_LS
    contrDu4e=error4eP1Energy(c4n,n4e,gradExact1,uApprox(:,1)) + ...
                        error4eP1Energy(c4n,n4e,gradExact2,uApprox(:,2));
                                        % L2-norm of Du - Du_LS
    
                                        
    %% Compute other local contributions to the LS functional
    for j=1:nrElems
        % local data
        nodes=n4e(j,:);         % nodes of this element
        coord=c4n(nodes,:)';    % coordinates of the nodes of this element
        sides=s4e(j,:);         % sides of this element
        area=area4e(j);         % area of this element

        signs=diag(normal4s(sides,:)*normal4e(:,:,j)')';
                                % scalar product of outer normals of this
                                % element and its sides

        % local coefficients
        sigmaloc=reshape(sigmaApprox(sides,:),6,1);
        
        % auxiliary variables
        N=coord(:)*ones(1,3)-repmat(coord(:,[3 1 2]),3,1);
        L=diag(signs.*[norm(N([3,4],2)),norm(N([5,6],3)),...
                                norm(N([1,2],1))]);

        % computation and assembly of the local stiffness matrices
        A1D1=(L*N'*(MA)*N*L)/(48*area);
        A=[A1D1,zeros(3);zeros(3),A1D1];
        
        % computation of L2-norm of sigma - sigma_LS
        intSigma1Baryc=[intSigma1Baryc1forE(j,:),...
                    intSigma1Baryc2forE(j,:),intSigma1Baryc3forE(j,:)];
        intSigma2Baryc=[intSigma2Baryc1forE(j,:),...
                    intSigma2Baryc2forE(j,:),intSigma2Baryc3forE(j,:)];
        
        H=[intSigma1Baryc*N(:,1),intSigma2Baryc*N(:,1);...
           intSigma1Baryc*N(:,2),intSigma2Baryc*N(:,2);...
           intSigma1Baryc*N(:,3),intSigma2Baryc*N(:,3)];
        
        prefactors=reshape(sigmaloc,3,2) .* repmat(diag(L),1,2) ./ (2*area);
        
        contrL2Sigma4e(j)=intSigmaSq4e(j) - 2*sum(diag(prefactors'*H))...
                                            + sigmaloc'*A*sigmaloc;
        
        % computation of L2-norm of div sigma - div sigma_LS
        Adiv1D=(L*ones(3)*L)/area;
        Adiv=[Adiv1D,zeros(3);zeros(3),Adiv1D];
        
        F=reshape(diag(L)*intF4e(j,:)./area,6,1);
        
        contrDivSigma4e(j)=sigmaloc'*Adiv*sigmaloc + ...
                            2*sigmaloc'*F + intFSq4e(j); 
        
        % local energy error
        error4e(j)=contrL2Sigma4e(j) + contrDivSigma4e(j) + contrDu4e(j);
    end
    
    err=sqrt(sum(error4e));
    
end


function [Lambda1,Lambda2,Lambda3]=barycentricCoords(c4n,n4e)
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
        % Barycentric coordinates and their gradients
        Lambda1=@(z)(Psi(:,1).*(z(:,1)-a3(:,1)) + Psi(:,3).*...
                                                    (z(:,2)-a3(:,2)));
        Lambda2=@(z)(Psi(:,2).*(z(:,1)-a3(:,1)) + Psi(:,4).*...
                                                    (z(:,2)-a3(:,2)));
        Lambda3=@(z)(1 - Lambda1(z) - Lambda2(z));
end
