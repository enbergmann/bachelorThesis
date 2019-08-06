function [err,errL2Sigma,errEnergy,error4e]=...
  error4eStokesLS(f,c4n,n4e,sigmaApprox,gradUApprox,sigmaExact,gradUExact,...
                  degreeF,degreeSigma,Adiv4e,F4e,area4e,contrDivSigma4e)
%% error4eStokesLS - compute error of LSFE solution for Stokes equations
% Computes the exact error of the LSFEM for the Stokes equations containing
% the H(div)-norm of the pseudostress and the energy norm of the velocity.
% This function benefits from parallel computing in the computation of the
% local error contributions.
%
% COMMAND:
%   [err,errL2Sigma,errEnergy,error4e]=...
%     error4eStokesLS(f,c4n,n4e,sigmaApprox,gradUApprox,sigmaExact,gradUExact,...
%                     degreeF,degreeSigma,Adiv4e,F4e,area4e,contrDivSigma4e)
%
% INPUT:
%   f       ... right-hand side of the problem definition
%   c4n     ... coordinates for the nodes of the mesh
%   n4e     ... nodes for the elements of the mesh
%   sigmaApprox basis coefficients of the numerical solution of
%               the pseudostress w.r.t. RT^2 basis
%   gradApprox. gradient of the discrete velocity
%   sigmaExact. function handle for the exact pseudostress
%   gradExact.. function handle for the gradient of the exact velocity
%   degreeF ... degree for the Gauss integration of f
%   degreeSigma degree for the Gauss integration of the sigmaExact
%   Adiv4e  ... local contributions to the stiffness matrix
%   F4e     ... local contributions to the RHS
%   area4e  ... area of each triangle (OPTIONAL)
%   contrDivSigma4e local contributions to the L2 norm of the divergence of
%               the pseudostress error (OPTIONAL)
%
% OUTPUT:
%   err     ... total error
%   errL2Sigma. error contribution to the L2 norm for each element
%   errEnergy.. error contribution to the energy norm for each element
%   error4e ... error for each element

% Copyright 2014 Philipp Bringmann
% Copyright 2011 Dietmar Gallistl (barycentricCoords subfunction)
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

%% PROCEED INPUT
if nargin<12; area4e=computeArea4e(c4n,n4e); end
if nargin<13
  computeDiv=true;
else
  computeDiv=false;
end

%% INITIALIZATION
% geometry
n4s=computeN4s(n4e);
s4e=computeS4e(n4e);
signs4e=computeSigns4e(c4n,n4e,n4s,s4e);
length4s=computeLength4s(c4n,n4s);
length4e=length4s(s4e);
coord4e=[c4n(n4e(:,1),:),c4n(n4e(:,2),:),c4n(n4e(:,3),:)];
% constants
nrElems=size(n4e,1);
% local coefficients
sigma4e=reshape(sigmaApprox(s4e,:),nrElems,6);

% functions
[Lambda1,Lambda2,Lambda3]=barycentricCoords(c4n,n4e);

% integrals
intSigmaBaryc1forE=...
    integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)sigmaExact(Gpts4p).*...
              repmat(Lambda1(Gpts4p),1,4),degreeSigma+1,area4e);
intSigmaBaryc2forE=...
    integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact(Gpts4p).*...
              repmat(Lambda2(Gpts4p),1,4)),degreeSigma+1,area4e);
intSigmaBaryc3forE=...
    integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)(sigmaExact(Gpts4p).*...
              repmat(Lambda3(Gpts4p),1,4)),degreeSigma+1,area4e);

% reordering variables to optimize parallelization
intSigma1Baryc1forE=intSigmaBaryc1forE(:,1:2);
intSigma1Baryc2forE=intSigmaBaryc2forE(:,1:2);
intSigma1Baryc3forE=intSigmaBaryc3forE(:,1:2);
intSigma2Baryc1forE=intSigmaBaryc1forE(:,3:4);
intSigma2Baryc2forE=intSigmaBaryc2forE(:,3:4);
intSigma2Baryc3forE=intSigmaBaryc3forE(:,3:4);

intSigmaSq4e=integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                      sum(sigmaExact(Gpts4p).*sigmaExact(Gpts4p),2),...
                      2*degreeSigma,area4e);

intFSq4e=integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)sum(f(Gpts4p).^2,2),...
                   2*degreeF,area4e);

%% CONSTANTS
MA=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],...
    [-4,-2,0,2,4],6,6);

%% VARIABLES
contrL2Sigma4e=zeros(nrElems,1);  % L2-norm of sigma - sigma_LS
if computeDiv
  contrDivSigma4e=zeros(nrElems,1); % L2-norm of div(sigma - sigma_LS)
end

%% COMPUTE ENERGY NORM OF VELOCITY ERROR
% contributions of the L2-norm of Du - Du_LS
contrDu4e=integrate(c4n,n4e,@(n4p,Gpts4p,Gpts4ref)...
                    sum((gradUExact(Gpts4p)-gradUApprox).^2,2),...
                    degreeSigma,area4e);

%% COMPUTE H(DIV) NORM OF PSEUDOSTRESS ERROR
% contributions of the H(div)-norm of sigma - sigma_LS
parfor j=1:nrElems
  %
  % PARALLIZABLE
  %

  % local data of each element
  coord=reshape(coord4e(j,:),2,3);
  area=area4e(j);

  % local coefficients
  sigmaloc=reshape(sigma4e(j,:),6,1);

  % auxiliary variables
  N=repmat(coord(:),1,3)-repmat(coord(:,[3 1 2]),3,1);
  L=diag(signs4e(j,:).*length4e(j,:));

  % computation and assembly of the local stiffness matrices
  A1D1=(L*N'*MA*N*L)/(48*area);
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
  if computeDiv
    contrDivSigma4e(j)=sigmaloc'*Adiv4e(:,:,j)*sigmaloc+2*sigmaloc'*F4e(:,j)...
                        +intFSq4e(j);
  end
end

%% ACCUMULATE
error4e=contrL2Sigma4e+contrDivSigma4e+contrDu4e;
err=sqrt(sum(error4e));
errL2Sigma=sqrt(sum(contrL2Sigma4e));
errEnergy=sqrt(sum(contrDu4e));

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
