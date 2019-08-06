function jump4s = P0NormalJump(c4n,n4e,n4sDb,n4s,s4e,sigma4e)
%% P0NormalJump - computes jumps in normal direction
%   c4n,n4e,n4sDb,n4sNb: mesh data
%   sigma4e: the values of a P0 function given for each element
%   g: a function giving values for the Neumann boundary
%
%   jump4s: the jumps of v across each side, for Neumann boundary
%           sides, g is the expected value.
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
% WARNING: modified, now without Neumann boundary

    %% Initialisation
    s4n = computeS4n(n4e);
    normal4e = computeNormal4e(c4n,n4e);
    jump4s = zeros(size(n4s,1),1);
    
    %% inner jumps
    for elem = 1 : size(n4e,1)
        jump4s(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + normal4e(:,:,elem)*sigma4e(elem,:)';
    end

    %% Dirichlet jumps = 0
    jump4s(diag(s4n(n4sDb(:,1),n4sDb(:,2)))) = 0;
    
    jump4s = abs(jump4s);

end
