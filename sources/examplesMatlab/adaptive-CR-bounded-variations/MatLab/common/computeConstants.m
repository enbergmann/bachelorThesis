function [cp,cF,kappa,capx2,c64,c7,c8,c9,c10,c11,c12,c13,c15,c16,c19]...
            = computeConstants(c1,c2,c3,c4,c5,c6,p,r,s,t,m,...
              Omega,diam,omega,Mpatch,M,sideCube,NormFPowerQ,c64)
%% computeConstants - compute constants to obtain guaranteed upper bounds
%
% Input:     c1,c2,c3,c4,c5,c6  constants related to the energy density     
%            p,r,s,t            parameters
%            m                  number of components of f
%            Omega,diam,omega   geometry of the domain and triangulation
%            Mpatch,M,sideCube
%            NormFPowerQ        Lq norm of f power q
% Optional:  c64                
%
% Output:    see below
          
    % conjugate exponent
    q = p/(p-1);
    % Poincare constant
    if p == 2
        cp = 1/pi;
    else
        cp = m^(1/2-1/p)*p*sin(pi/p)/(pi*(p-1)^(1/p));
    end
    % Friedrichs constant
    cF = sideCube;
    % nonconforming interpolation
    kappa = m^(1/2-1/p)*cp + 2^(1/q)/(sqrt(3)*2*(q+2)^(1/q));
    % enrichment operator
    capx = 3^(p+1/2)*cot(omega)*max([Mpatch^(1-p/2), Mpatch^(p/2 - 1)]);
    capx = capx/(2^((p-2)/2)*(2+p/(p-1))^(p-1)...
            *sin(omega)^(M*max([p-2,2-p]))*(1-cos(pi/Mpatch))^(p/2));
    capx = capx^(1/p);
    % conforming companion of n-th order
    capx2 = computeFactorApprox(p)*capx;
    % c7
    c7 = (q/c1*(c2+c4)*Omega + cF^q/c1^q*NormFPowerQ)^(1/p);
    % c8
    c8 = 1/((c3*p)^(q/p)*q);
    % c9
    c9 = 1/((c1*p)^(q/p)*q);
    % c10
    c10 = (2^(q-1)*q*c9)^(p-1)*p*max([c2 + c4, 2^(q-1)*c9 - c8]);
    % if c64 is not given
    if nargin < 19
        c64 = c6*r*(1+max([1,3^(s-1)])*c10^(s/p)...
                *max([(1+3^((q-1)/p)*(2*c5)^(q/p))^s, 3^((q-1)*s/p)]));
    end
    % c11
    c11 = c8^(-1/q)*((c2 + c4)*Omega + cF^q/((c1*p)^(q/p)*q)*NormFPowerQ)^(1/q);
    % c12
    c12 = 2^(1/p)*c5*(Omega + c7^p)^(1/q);
    % c13
    c13 = c12 + sqrt(m)*diam*NormFPowerQ/3;
    % c15
    c15 = max([3,3^(t-1)])*c6*(Omega + 2*c7^p)^(t-1);
    % c16
    c16 = max([3,3^(t-1)])*c64*(Omega + 2*c7^p)^(t-1);
    % c19
    c19 = max([3,3^(t-1)])*c6*r*(Omega + c7^p + c10*(Omega + 2^(q-1)*c11^q + 2^(q-1)*c13^q))^(t-1);
end

function val =  computeFactorApprox(p)
    val = 0;
    m = floor(p)-1;
    for j = 0:m
        for k = 0:m-j
            val = val + factorial(m-k+1)*factorial(m-j)*factorial(j+k+1)/...
                    (factorial(j)*factorial(k)*factorial(m-j-k));
        end
    end
    val = 3^(m+1)*2^(p+1)*factorial(m)*val/factorial(4+2*m);
    val = 1 + 3^(1/2)*val^(1/p);
end

% Copyright 2018 Tran Ngoc Tien
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. No WARRANTY!