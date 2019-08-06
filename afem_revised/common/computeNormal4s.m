function normal4s = computeNormal4s(c4n,n4s)
%% computeNormal4s - Normals for sides.
%   computeNormal4s(c4n, n4s) computes the outer unit normal vectors for each
%                         side of the decomposition. For inner sides only
%                         one of the two possible normals is computed.
%                         c4n and n4s are as specified in the
%                         documentation.
%
%   See also: computeN4s, computeNormal4e, computeTangent4s
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
    
    if isempty(n4s)
        normal4s = zeros(0,2);
        return;
    end

    %% Compute length4s.
    c4start = c4n(n4s(:,1),:);
    c4end = c4n(n4s(:,2),:);
    length4s = sqrt(sum((c4end-c4start).^2,2));
    
    %% Compute tangent4s.              
    tangent4s = (c4end-c4start)./[length4s length4s];
    
    %% Compute normal4s.
    normal4s = [tangent4s(:,2), -tangent4s(:,1)];
end