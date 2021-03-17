function mid4e = computeMid4e(c4n, n4e)
%% computeMid4e - midpoints for elements.
%   computeMid4e(c4n, n4e) computes the midpoint for each element of the
%                     decomposition. c4n and n4e are as specified in the
%                     documentation.
%
%   See also: computeArea4e, computeMid4s
%
% This file is taken from the AFEM software package.
% (C) 2009 Num. Analysis group, Prof. Carstensen, HU Berlin
% Licensed under GNU GPL 3+. No warranty! See LICENSE.txt.

    if isempty(n4e)
        mid4e = zeros(0,2);
        return;
    end

    %% Compute mp4e.
    mid4e = ( c4n(n4e(:,1),:) + c4n(n4e(:,2),:) + c4n(n4e(:,3),:) ) / 3;
end