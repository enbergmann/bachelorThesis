function mid4s = computeMid4s(c4n, n4s)
%% computeMid4s - midpoints for sides.
%   computeMid4s(c4n, n4s) computes the midpoint for each side of the
%                     decomposition. c4n and n4s are as specified in the
%                     documentation.
%
%   See also: computeN4s, computeMid4e
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.
    
    if isempty(n4s)
        mid4s = zeros(0,2);
        return;
    end

    %% Compute mid4s.
    mid4s = 0.5 * ( c4n(n4s(:,1),:) + c4n(n4s(:,2),:) );
end