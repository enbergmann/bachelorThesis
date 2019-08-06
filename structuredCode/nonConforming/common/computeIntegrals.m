function [rhsInt1,rhsInt2,rhsInt3] = computeIntegrals(f,c4n,n4e,degree,area4e)
  % phi1 = @(x)( x(:,1) );
  % phi2 = @(x)( x(:,2) );
  % phi3 = @(x)( 1 - x(:,1) - x(:,2) );
  % 
  % psi1 = @(x)( 1-2*phi3(x) );
  % psi2 = @(x)( 1-2*phi1(x) );
  % psi3 = @(x)( 1-2*phi2(x) );
  
  %INT1 = @(n4p,Gpts4p,Gpts4ref) ...
  %    psi1(Gpts4ref).*f(Gpts4p);
  %INT2 = @(n4p,Gpts4p,Gpts4ref) ...
  %    psi2(Gpts4ref).*f(Gpts4p);        
  %INT3 = @(n4p,Gpts4p,Gpts4ref) ...
  %    psi3(Gpts4ref).*f(Gpts4p);
  %
  %temp1 = integrate(INT1, c4n, n4e, degree, area4e);
  %temp2 = integrate(INT2, c4n, n4e, degree, area4e);
  %temp3 = integrate(INT3, c4n, n4e, degree, area4e);

  [Lambda1,Lambda2,Lambda3] = barycentricCoords(c4n,n4e);
  rhsInt1 = ...
      integrate(@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p).*(1 - 2*Lambda3(Gpts4p)))...
    ,c4n,n4e,degree+1,area4e);
  rhsInt2 = ...
      integrate(@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p).*(1 - 2*Lambda1(Gpts4p)))...
    ,c4n,n4e,degree+1,area4e);
  rhsInt3 = ...
      integrate(@(n4p,Gpts4p,Gpts4ref)(f(Gpts4p).*(1 - 2*Lambda2(Gpts4p)))...
    ,c4n,n4e,degree+1,area4e);
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
