function uNC = interpolationNC(u,c4n,n4e,n4s)
  % computes the CR coefficients of a given function handle

  %% naive
  % mid4s = computeMid4s(c4n,n4s);
  % uNC = u(mid4s);

  %% actual CR interpolation 

  length4s = computeLength4s(c4n,n4s);

  INT = @(n4p,Gpts4p,Gpts4ref) ...
      u(Gpts4p);
  
  uNC = integrate(INT, c4n, n4s, 20, length4s)./length4s;
end
