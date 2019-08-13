function uNC = interpolationNC(u, c4n, n4e, n4s)
% Computes the Crouzeix-Raviart interpolation uNC (vanishing in the midpoints
% of boundary edges) with respect to the triangulation given by [c4n, n4e, n4s]
% of the function u.
%
% interpolationNC.m
% input:  u   - 'function handle' of the function to be interpolated
%         c4n - coordinates for nodes
%         n4e - nodes for elements
%         n4s - nodes for sides
%
% output: uNC - 'function handle' of the interpolation of u

  length4s = computeLength4s(c4n,n4s);

  integrand = @(n4p,Gpts4p,Gpts4ref) ...
      u(Gpts4p);
  
  uNC = integrate(integrand, c4n, n4s, 20, length4s)./length4s;
end







%% naive
% mid4s = computeMid4s(c4n,n4s);
% uNC = u(mid4s);
