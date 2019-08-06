%function [ contrDu4e, contrSigma4e, error4e, contrDivSigma4e, contrL2Sigma4e, contrL2U4e, contrU4e] = ...
%error4eElasticityLS(geometryData, u, sigma, problemData, lambda,mu ,exSolutionData, localData,solutionData, intDegree, method)
function [ contrDu4e, contrSigma4e, error4e, contrDivSigma4e, contrL2Sigma4e, contrL2U4e, contrU4e] = ...
    error4eElasticityLS(geometryData, problemData, exSolutionData, localData, solutionData, parameter);

u     = solutionData.u;
sigma = solutionData.sigma;

lambda    = parameter.lambda;
lame_mu   = parameter.lame_mu;
intDegree = parameter.intDegree;
method    = parameter.method;

% geometry
c4n = geometryData.c4n; %geometryData{1};
n4e = geometryData.n4e; %geometryData{2};
n4s = computeN4s(n4e);

nrElems = size(n4e,1);
nrNodes = size(c4n,1);
nrSides  = size(n4s,1);

f = problemData.f; %problemData{1};

u1Cex        = exSolutionData.u1Cex; %exSolutionData{1};
u2Cex        = exSolutionData.u2Cex; %exSolutionData{2};
gradUexact1C = exSolutionData.gradUexact1C; %exSolutionData{8};
gradUexact2C = exSolutionData.gradUexact2C; %exSolutionData{9};
sigmaExact   = exSolutionData.sigmaExact; %exSolutionData{7};

% read local Matrizies
L4e     = localData.L4e; %localData{3};
N4e     = localData.N4e; %localData{4};
F4e     = localData.F4e; %localData{8};
Adiv4e  = localData.Adiv4e; %localData{9};
AL24e   = localData.AL24e; %localData{11};
area4e  = localData.area4e; %localData{15};
sigma4e = solutionData.sigma4e; %solutionData{10};

%% integrals
% integral of each component of f over each element
  %intF4e = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p)),deg,area4e);
% integral of the scalar product f^2 over each element
%intFSq4e     = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref) sum(f(Gpts4p).^2,2),2*intDegree,[nrElems,1],area4e);
intFSq4e     = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref) sum(f(Gpts4p).^2,2),2*intDegree,area4e);
% integral of sigmaExact : sigmaExact
intSigmaSq4e = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                            permute(sum(sum(sigmaExact(Gpts4p).^2)),[3 1 2]),...
                            intDegree,[nrElems,1,1],area4e);
                           %sum(sum(SigmaExact(Gpts4p,lambda,mu).^2))
% intSigmaSq4e = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
%                             permute(sum(sum(sigmaExact(Gpts4p).^2)),[3 1 2]),...
%                             intDegree,area4e);
%                            %sum(sum(SigmaExact(Gpts4p,lambda,mu).^2))                           


[Lambda1,Lambda2,Lambda3]=barycentricCoords(c4n,n4e);

intSigmaBaryc1forE = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                     (reshape(permute(sigmaExact(Gpts4p),[3,2,1]),size(n4p,1),4,1).*repmat(Lambda1(Gpts4p),1,4)),...
                     intDegree+1,[nrElems,4,1],area4e);
intSigmaBaryc2forE = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                     (reshape(permute(sigmaExact(Gpts4p),[3,2,1]),size(n4p,1),4,1).*repmat(Lambda2(Gpts4p),1,4)),...
                     intDegree+1,[nrElems,4,1],area4e);
intSigmaBaryc3forE = parIntegrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                     (reshape(permute(sigmaExact(Gpts4p),[3,2,1]),size(n4p,1),4,1).*repmat(Lambda3(Gpts4p),1,4)),...
                     intDegree+1,[nrElems,4,1],area4e);
% intSigmaBaryc1forE = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
%                      (reshape(permute(sigmaExact(Gpts4p),[3,2,1]),size(n4p,1),4,1).*repmat(Lambda1(Gpts4p),1,4)),...
%                      intDegree+1,area4e);
% intSigmaBaryc2forE = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
%                      (reshape(permute(sigmaExact(Gpts4p),[3,2,1]),size(n4p,1),4,1).*repmat(Lambda2(Gpts4p),1,4)),...
%                      intDegree+1,area4e);
% intSigmaBaryc3forE = integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
%                      (reshape(permute(sigmaExact(Gpts4p),[3,2,1]),size(n4p,1),4,1).*repmat(Lambda3(Gpts4p),1,4)),...
%                      intDegree+1,area4e);

%% reordering variables to optimize parallelization
intSigma1Baryc1forE = intSigmaBaryc1forE(:,1:2);
intSigma1Baryc2forE = intSigmaBaryc2forE(:,1:2);
intSigma1Baryc3forE = intSigmaBaryc3forE(:,1:2);

intSigma2Baryc1forE = intSigmaBaryc1forE(:,3:4);
intSigma2Baryc2forE = intSigmaBaryc2forE(:,3:4);
intSigma2Baryc3forE = intSigmaBaryc3forE(:,3:4);

% preallocation
contrDivSigma4e = zeros(nrElems,1);
contrL2Sigma4e  = zeros(nrElems,1);
contrSigma4e    = zeros(nrElems,1);
error4e         = zeros(nrElems,1);

switch method
    case 'S1'
        contrL2U4e = error4eP1L2(c4n,n4e,u1Cex, u(1:2:2*nrNodes)) + ...
                     error4eP1L2(c4n,n4e,u2Cex, u(2:2:2*nrNodes));
        
        contrDu4e  = error4eP1Energy(c4n,n4e,gradUexact1C,u(1:2:2*nrNodes)) + ...
                     error4eP1Energy(c4n,n4e,gradUexact2C,u(2:2:2*nrNodes));
    case 'KS'
        contrL2U4e = error4eP1L2(c4n,n4e,u1Cex, u(1:nrNodes)) + ...
                     error4eCRL2(c4n,n4e,u2Cex, u(nrNodes+1:nrNodes+nrSides));

        contrDu4e  = error4eP1Energy(c4n,n4e,gradUexact1C,u(1:nrNodes)) + ...
                     error4eCREnergy(c4n,n4e,gradUexact2C,u(nrNodes+1:nrNodes + nrSides));
    case 'CR'
        contrL2U4e = error4eCRL2(c4n,n4e,u1Cex, u(1:2:2*nrSides)) + ...
                     error4eCRL2(c4n,n4e,u2Cex, u(2:2:2*nrSides));

        contrDu4e  = error4eCREnergy(c4n,n4e,gradUexact1C,u(1:2:2*nrSides)) + ...
                     error4eCREnergy(c4n,n4e,gradUexact2C,u(2:2:2*nrSides));
end

if matlabpool('size') == 0
  talk = 1;
else
  talk = 0;
end

%% LOOP
parfor j = 1: nrElems
    if talk
      fprintf('\r                                        ')%replaces current line by 40 spaces
      fprintf('\r ... elements to compute %d/%d',j,nrElems)
    end
    % H(div) norm of sigma
    % divergence contributiuon
    contrDivSigma4e(j) = sigma4e(j,:)*Adiv4e(:,:,j)*sigma4e(j,:)'- ...     
                         2.*sigma4e(j,:)*F4e(:,j) + ...
                         intFSq4e(j);
                         
    % L2 contribution
    intSigma1Baryc=[intSigma1Baryc1forE(j,:),...
              intSigma1Baryc2forE(j,:),intSigma1Baryc3forE(j,:)];
    intSigma2Baryc=[intSigma2Baryc1forE(j,:),...
              intSigma2Baryc2forE(j,:),intSigma2Baryc3forE(j,:)];
    
    N    = N4e(:,:,j);
    L    = L4e(:,:,j);
    area = area4e(j);
    
    H=[intSigma1Baryc*N(:,1),intSigma2Baryc*N(:,1);...
       intSigma1Baryc*N(:,2),intSigma2Baryc*N(:,2);...
       intSigma1Baryc*N(:,3),intSigma2Baryc*N(:,3)];

    
    %sigmaloc=reshape(sigma4e(j,:),6,1);
    %prefactors=reshape(sigmaloc,3,2) .* repmat(diag(L),1,2) ./ (2*area);
    prefactors=reshape(sigma4e(j,:)',2,3)' .* repmat(diag(L),1,2) ./ (2*area);
    
    contrL2Sigma4e(j) = intSigmaSq4e(j)            - ...                               
                        2*sum(diag(prefactors'*H)) + ...        
                        sigma4e(j,:)*AL24e(:,:,j)*sigma4e(j,:)';
end

if talk
  fprintf('\r                                        \r')%replaces current line by 40 spaces
end

contrU4e     = contrL2U4e      + contrDu4e;
contrSigma4e = contrDivSigma4e + contrL2Sigma4e;
error4e      = contrSigma4e    + contrU4e;

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
% Barycentric coordinates and their gradients
Lambda1=@(z)(Psi(:,1).*(z(:,1)-a3(:,1)) + Psi(:,3).*...
                                            (z(:,2)-a3(:,2)));
Lambda2=@(z)(Psi(:,2).*(z(:,1)-a3(:,1)) + Psi(:,4).*...
                                            (z(:,2)-a3(:,2)));
Lambda3=@(z)(1 - Lambda1(z) - Lambda2(z));
end

