function osc4e = computeOscError(c4n, n4e, f, m, q, area4e, intF4e)
%% computeOscError - compute oscillation
%
% Input:     c4n            coordinates for the nodes of the mesh
%            n4e            nodes for the elements of the mesh
%            f              body force
%            m              number of components of f
%            q              parameter
%            intF4e         integral of F for each element
%            area4e         area of each element
%
% Output:    osc4e          oscillation of f

    %% COMPUTATION
    % Compute integral mean
    mean4e = intF4e./repmat(area4e,1,m);
    
    % bestapproximation error
    osc4e = repmat(area4e,1,m).^(q/2).*integrate(c4n,n4e,...
                 @(n4p,Gpts4p,Gpts4ref)((sum((f(Gpts4p) - mean4e).^2,2)).^(q/2)),10,area4e);
end

% Copyright 2018 Tran Ngoc Tien
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. No WARRANTY!