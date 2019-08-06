function length4s = computeLength4s(c4n,n4s)
%% computeLength4s - Lengths for sides.
%   computeLength4s(c4n, n4s) computes the length of each side of the
%                         decomposition. c4n and n4s are as specified in the
%                         documentation.
%
%   See also: computeN4s, computeArea4e
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
    
    if isempty(n4s)
        length4s = zeros(0,1);
        return;
    end

    %% Compute length4s in a vectorised manner.
    length4s = sqrt( sum( (c4n(n4s(:,2),:) - c4n(n4s(:,1),:)).^2, 2) );
end