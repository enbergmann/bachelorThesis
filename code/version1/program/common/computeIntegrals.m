function [temp1,temp2,temp3] = computeIntegrals(f,c4n,n4e,degree,area4e)
  phi1 = @(x)( x(:,1) );
  phi2 = @(x)( x(:,2) );
  phi3 = @(x)( 1 - x(:,1) - x(:,2) );
  
  psi1 = @(x)( 1-2*phi3(x) );
  psi2 = @(x)( 1-2*phi1(x) );
  psi3 = @(x)( 1-2*phi2(x) );
  
  INT1 = @(n4p,Gpts4p,Gpts4ref) ...
      psi1(Gpts4ref).*f(Gpts4p);
  INT2 = @(n4p,Gpts4p,Gpts4ref) ...
      psi2(Gpts4ref).*f(Gpts4p);        
  INT3 = @(n4p,Gpts4p,Gpts4ref) ...
      psi3(Gpts4ref).*f(Gpts4p);
  
  temp1 = integrate(INT1, c4n, n4e, degree, area4e);
  temp2 = integrate(INT2, c4n, n4e, degree, area4e);
  temp3 = integrate(INT3, c4n, n4e, degree, area4e);
end

