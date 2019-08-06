function [wNormF4e, wNormP0F4e]...
            = computeWeightedNormF(c4n, n4e, f, m, q, r, t, intF4e, degreeF, area4e)
%% computeWeightedNormF - compute norm of hf
%
% Input:     c4n            coordinates for the nodes of the mesh
%            n4e            nodes for the elements of the mesh
%            f              body force
%            m              number of components of f
%            q,r,t          parameters
%            intF4e         integral of F for each element
%            degreeF        degree for Gauss integration of f
%            area4e         area of each element
%
% Output:    wNormF4e       Lq norm of hf
%            wNormP0F4e     L r/t norm of h Pi0 f

    degree=3*degreeF;
    wNormF4e = area4e.^(q/2).*integrate(c4n,n4e,...
                 @(n4p,Gpts4p,Gpts4ref)((sum(f(Gpts4p).^2,2)).^(q/2)),degree);
    P0F4e = intF4e./repmat(area4e,1,m);
    wNormP0F4e = area4e.^(1 + r/(2*t)).*sum(P0F4e.^2,2).^(r/(2*t));
end

% Copyright 2018 Tran Ngoc Tien
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. No WARRANTY!