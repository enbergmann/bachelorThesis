function error4e = error4eP1L2(c4n,n4e,uExact,uApprox)
%% error4eP1L2 - compute the exact error on elements
%  Input:  c4n,n4e     - mesh
%          uExact      - the exact function
%          uApprox     - AFEM P1 solution
%
%  Output: error4e     - the exact squared error on each element

%% Compute error
  error4e = (integrate(c4n,n4e,@(n4p, Gpts4p, Gpts4ref) (...
           uExact(Gpts4p)-...
            ((1 - Gpts4ref(:,1) - Gpts4ref(:,2))*uApprox(n4p(:,1)) +...
            (Gpts4ref(:,1))*uApprox(n4p(:,2)) +...
            (Gpts4ref(:,2))*uApprox(n4p(:,3)))).^2,4));      
        
end
